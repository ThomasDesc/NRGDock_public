import os

import numpy as np
from numba import njit
from os import walk
from os import path
from complementarity_function import get_cf
import math as m
import datetime
import sys


def njit(f):
    return f


def load_rad_list(filepath):
    dict = {}
    with open(filepath) as f:
        lines = f.readlines()
        for a, line in enumerate(lines):
            line = line.strip().split(' ')
            dict[line[0]] = [a+1, float(line[1])]
    return dict


@njit
def get_radius_number(letter_type, rad_dict):
    if isinstance(letter_type, int):
        return list(rad_dict.values())[letter_type-1][1]
    else:
        letter_type = letter_type.upper().replace(' ', '')
        return rad_dict[letter_type][0], rad_dict[letter_type][1]



@njit
def get_params_dict(config_file):
    """ Params will be a dictionary of the important parameters, loaded from a text file. For example, number of
    orientations per dimension for the ligand, number of points in the binding site grid, etc."""
    params_dict = dict()
    with open(config_file) as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('N_ORIENTATIONS'):
            params_dict['N_ORIENTATIONS'] = int(line.split()[-1])
        elif line.startswith('WATER_RADIUS'):
            params_dict['WATER_RADIUS'] = float(line.split()[-1])
        elif line.startswith('DOT_DIVISION'):
            params_dict['DOT_DIVISION'] = float(line.split()[-1])
        elif line.startswith('GRID_PLACEHOLDER'):
            params_dict['GRID_PLACEHOLDER'] = float(line.split()[-1])
        elif line.startswith('DEVICE'):
            params_dict['DEVICE'] = line.split()[-1]
        elif line.startswith('GA_SPHERE_RADIUS'):
            params_dict["GA_SPHERE_RADIUS"] = float(line.split()[-1])
        elif line.startswith('GA_GENERATIONS'):
            params_dict["GA_GENERATIONS"] = int(line.split()[-1])
        elif line.startswith('OUTPUT_PDB'):
            params_dict["OUTPUT_PDB"] = line.split()[-1]
    return params_dict


@njit
def build_ligand_list(ligand_path):
    ligand_list = []
    _, _, filenames = next(walk(ligand_path), (None, None, []))
    for filename in filenames:
        if filename.endswith(".mol2"):
            ligand_list.append(ligand_path + filename)
    return ligand_list


@njit
def load_atoms_mol2(filename, ligand_line_start, ligand_last_line, rad_dict):
    coord_start = 0
    with open(filename) as f:
        lines = f.readlines()
    if ligand_line_start and ligand_last_line is not None:
        lines = lines[int(ligand_line_start):int(ligand_last_line)]
    n_atoms = 0
    for line in lines:
        if line.startswith('@<TRIPOS>ATOM'):
            coord_start = 1
        if line.startswith('@<TRIPOS>BOND') or line.startswith("@<TRIPOS>UNITY"):
            break
        if coord_start == 1 and line.startswith('@<TRIPOS>ATOM') == False and line[8] != 'H':
            n_atoms += 1
    atoms_xyz = np.zeros((n_atoms, 3), dtype=np.float32)
    atoms_numbers_types = np.zeros((n_atoms, 2), dtype=np.float32)  # 1st column is number, 2nd will be atom type
    atoms_radius = np.zeros((n_atoms, 1), dtype=np.float32)
    counter = 0
    coord_start = 0
    for line in lines:
        if line.startswith('@<TRIPOS>ATOM'):
            coord_start = 1
        if line.startswith('@<TRIPOS>BOND') or line.startswith("@<TRIPOS>UNITY"):
            break
        if coord_start == 1 and line.startswith('@<TRIPOS>ATOM') == False and line[8] != 'H':
            atoms_xyz[counter] = np.array([float(line[-60:-52]), float(line[-50:-42]), float(line[-40:-32])])
            atoms_numbers_types[counter][0] = int(line[:7])
            atoms_numbers_types[counter][1], atoms_radius[counter] = get_radius_number(line[46:52], rad_dict)
            counter += 1
    atoms_numbers_types_sorted = atoms_numbers_types[atoms_numbers_types[:, 0].argsort()]
    return atoms_xyz, atoms_numbers_types_sorted, atoms_radius


@njit
def load_binding_site_grid(params, binding_site):
    """ This will load the binding site spheres, find the extreme values for x, y, z, then build a grid of points
    and remove the ones that clash with the target atoms.
    """
    dot_division = params['DOT_DIVISION']
    list = []
    box_limit = []
    grid_coords = []
    with open(binding_site) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('ATOM'):
                temp_array = []
                line = line.strip()
                line = line.split(" ")
                for item in line:
                    if len(temp_array) > 2:
                        break
                    else:
                        if len(item.split('.')) > 1:
                            temp_array.append(float(item))
                temp_array.append(float(line[-1]))     # sphere radius
                list.append(temp_array)
    a = np.array(list)
    x = [round(np.min(a[:, 0]) - a[np.argmin(a[:, 0]), 3], 3), round(np.max(a[:, 0]) + a[np.argmax(a[:, 0]), 3], 3)]  # find min and max coords for box as well as remove or add sphere radius
    y = [round(np.min(a[:, 1]) - a[np.argmin(a[:, 1]), 3], 3), round(np.max(a[:, 1]) + a[np.argmax(a[:, 1]), 3], 3)]
    z = [round(np.min(a[:, 2]) - a[np.argmin(a[:, 2]), 3], 3), round(np.max(a[:, 2]) + a[np.argmax(a[:, 2]), 3], 3)]

    for dot_x in np.arange(x[0], x[1], dot_division):
        for dot_y in np.arange(y[0], y[1], dot_division):
            for dot_z in np.arange(z[0], z[1], dot_division):
                coords = np.array([round(dot_x, 3), round(dot_y, 3), round(dot_z, 3)])
                for row in a:
                    distance = np.linalg.norm(coords - row[:3])
                    if distance < row[3]:
                        grid_coords.append(coords)
                        break
    #write_test(grid_coords, "OG_grid", PATH)
    return grid_coords


@njit
def write_test(coord_list, name, path):
    textfile = open(f"{path}{name}.pdb", 'w')
    counter = 7000
    for line in coord_list:
        textfile.write("ATOM   {:>4}  C   DOT X   1{:>12} {:>7} {:>7}   1.00  0.10 \n".format(str(counter),
                                                                                              str(round(line[0], 3)),
                                                                                              str(round(line[1], 3)),
                                                                                              str(round(line[2], 3))))
        counter += 1
    textfile.close()


@njit
def build_3d_cube_grid(params, target_atoms_xyz, atoms_radius, cw_factor=4):
    """ This will scan all the target atoms and determine the extreme coordinates of the grid where all atoms will be
    indexed, and the length of the cubes. It will also index the target atoms in a grid using these params"""

    water_radius = params['WATER_RADIUS']
    grid_placeholder = params['GRID_PLACEHOLDER']

    max_rad = np.amax(atoms_radius, axis=0)
    cell_width = 2 * (max_rad + water_radius)

    # max_xyz = np.amax(target_atoms_xyz, axis=0)
    max_xyz = np.zeros(3)
    max_xyz[0] = np.max(target_atoms_xyz[:, 0]) + cell_width*cw_factor
    max_xyz[1] = np.max(target_atoms_xyz[:, 1]) + cell_width*cw_factor
    max_xyz[2] = np.max(target_atoms_xyz[:, 2]) + cell_width*cw_factor

    # min_xyz = np.amin(target_atoms_xyz, axis=0)
    min_xyz = np.zeros(3)
    min_xyz[0] = np.min(target_atoms_xyz[:, 0]) - cell_width*cw_factor
    min_xyz[1] = np.min(target_atoms_xyz[:, 1]) - cell_width*cw_factor
    min_xyz[2] = np.min(target_atoms_xyz[:, 2]) - cell_width*cw_factor

    lengths = ((max_xyz - min_xyz) / cell_width).astype(np.int32) + 1
    temp_grid = []
    for i in range(lengths[0]):
        temp_grid.append([])
        for j in range(lengths[1]):
            temp_grid[i].append([])
            for k in range(lengths[2]):
                temp_grid[i][j].append([])
    for i, row in enumerate(target_atoms_xyz):
        grid_indices = ((row[:3] - min_xyz) / cell_width).astype(np.int32)
        temp_grid[grid_indices[0]][grid_indices[1]][grid_indices[2]].append(i)
    max_cell_len = 0
    for row in temp_grid:
        for col in row:
            for cell in col:
                n = len(cell)
                if n > max_cell_len:
                    max_cell_len = n
    grid = np.full((lengths[0], lengths[1], lengths[2], max_cell_len), grid_placeholder, dtype=np.int32)
    for i in range(lengths[0]):
        for j in range(lengths[1]):
            for k in range(lengths[2]):
                for x, v in enumerate(temp_grid[i][j][k]):
                    grid[i][j][k][x] = v
    return grid, min_xyz, cell_width



@njit
def load_energy_matrix(filename):
    with open(filename) as f:
        lines = f.readlines()
    energy_matrix = np.zeros((41, 41), dtype=np.float32)
    for line in lines:
        line = line.replace(" ","")
        line = line.strip()
        line = line.split('=')
        line_arr = line[0].split('-')
        line_arr.append(line[1])
        energy_matrix[int(line_arr[0]), int(line_arr[1])] = line_arr[2]
        energy_matrix[int(line_arr[1]), int(line_arr[0])] = line_arr[2]
        pass

    return energy_matrix


@njit
def center_coords(ligand_atoms_xyz):
    length = ligand_atoms_xyz.shape[0]
    sum_x = np.sum(ligand_atoms_xyz[:, 0])/length
    sum_y = np.sum(ligand_atoms_xyz[:, 1])/length
    sum_z = np.sum(ligand_atoms_xyz[:, 2])/length

    centroid_coords = np.array([sum_x, sum_y, sum_z])
    for i in range(len(ligand_atoms_xyz)):
        ligand_atoms_xyz[i] -= centroid_coords

    return ligand_atoms_xyz

@njit
def Rx(theta):
    return np.matrix([[1, 0, 0],
                      [0, m.cos(theta), -m.sin(theta)],
                      [0, m.sin(theta), m.cos(theta)]])

@njit
def Ry(theta):
    return np.matrix([[m.cos(theta), 0, m.sin(theta)],
                      [0, 1, 0],
                      [-m.sin(theta), 0, m.cos(theta)]])

@njit
def Rz(theta):
    return np.matrix([[m.cos(theta), -m.sin(theta), 0],
                      [m.sin(theta), m.cos(theta), 0],
                      [0, 0, 1]])


@njit
def rotate_ligand(ligand_atoms_xyz, params):
    ligand_atoms_xyz = center_coords(ligand_atoms_xyz)
    divisions = params['N_ORIENTATIONS']
    single_rotation = 360/divisions
    rotation_counter = 0
    ligand_orientations = np.zeros((divisions**3, len(ligand_atoms_xyz), 3), dtype=np.float32)
    for x in range(divisions):
        for y in range(divisions):
            for z in range(divisions):
                x = m.radians(x*single_rotation)
                y = m.radians(y*single_rotation)
                z = m.radians(z*single_rotation)
                rotation_matrix = Rx(x) * Ry(y) * Rz(z)
                for i, coord in enumerate(ligand_atoms_xyz):
                    ligand_orientations[rotation_counter][i] = np.array(np.dot(rotation_matrix, coord))
                #write_test(ligand_orientations[rotation_counter], "ligand_rotation_" + str(rotation_counter))
                rotation_counter += 1
    return ligand_orientations


def import_pred_list(name):
    with open("./predictor_list/" + name) as f:
        lines = f.readlines()
        surface_array = np.zeros(len(lines), dtype=np.float64)
        for i, line in enumerate(lines):
            surface_array[i] = float(line)
        return surface_array


def clean_bindig_site_grid(params, target_grid, binding_site_grid, min_xyz, cell_width, target_atoms_xyz, verbose):
    index = []
    for a, point in enumerate(binding_site_grid):
        grid_index = ((point - min_xyz) / cell_width).astype(np.int32)
        for i_offset in [-1, 0, 1]:
            for j_offset in [-1, 0, 1]:
                for k_offset in [-1, 0, 1]:
                    i = i_offset + grid_index[0]
                    j = j_offset + grid_index[1]
                    k = k_offset + grid_index[2]
                    for neighbour in target_grid[i][j][k]:
                        if neighbour == -1:
                            break
                        else:
                            dist = np.sqrt((target_atoms_xyz[neighbour][0] - point[0]) ** 2 + (target_atoms_xyz[neighbour][1] - point[1]) ** 2 + (target_atoms_xyz[neighbour][2] - point[2]) ** 2)
                            if dist <= 2.0:
                                index.append(a)
                                break
    cleaned_binding_site_grid = np.delete(binding_site_grid, index, 0)
    if verbose:
        print("REMARK Deleted binding site grid dots: ", len(index))
        print("REMARK Total binding site grid dots: ", len(cleaned_binding_site_grid))
        print("REMARK Total CF evaluations per ligand: ", len(cleaned_binding_site_grid) * params["N_ORIENTATIONS"]**3)
    #write_test(cleaned_binding_site_grid, "cleaned", PATH)
    return cleaned_binding_site_grid


@njit
def main(config_file, ligands_list, binding_site, target, energy_matrix_file, last_line):
    params_dict = get_params_dict(config_file)
    print("REMARK software: main")
    print(f"REMARK orientations: {params_dict['N_ORIENTATIONS']}")
    print(f"REMARK dot separation: {params_dict['DOT_DIVISION']}")
    pred_05 = import_pred_list("pred_05.txt")
    pred_95 = import_pred_list("pred_95.txt")
    rad_dict = load_rad_list("./radius_list.txt")
    #target_atoms_xyz, _, target_atoms_numbers_types, atoms_radius = load_atoms(target, rad_dict)
    target_atoms_xyz, target_atoms_numbers_types, atoms_radius = load_atoms_mol2(target, None, None, rad_dict)
    target_grid, min_xyz, cell_width = build_3d_cube_grid(params_dict, target_atoms_xyz, atoms_radius)
    original_grid = load_binding_site_grid(params_dict, binding_site)
    binding_site_grid = clean_bindig_site_grid(params_dict, target_grid, original_grid, min_xyz, cell_width, target_atoms_xyz, verbose)
    energy_matrix = load_energy_matrix(energy_matrix_file)
    n_cf_evals = len(binding_site_grid) * params_dict["N_ORIENTATIONS"]**3
    cfs_list_by_ligand = np.zeros((len(ligands_list)), dtype=np.float32)
    non_zero_list = []
    non_clash_list = []
    total_evals_list = []
    if verbose and Time:
        print(f"REMARK ligand start: {datetime.datetime.now().time()}")
    for i, ligand in enumerate(ligands_list):
        cfs_list = np.zeros((n_cf_evals, 3), dtype=np.float32)
        ligand_atoms_xyz, ligand_atoms_numbers_types, ligand_atoms_radius = load_atoms_mol2(path_to_ligands, ligand, last_line, rad_dict)
        atm_quantity.append(len(ligand_atoms_xyz))
        ligand_orientations = rotate_ligand(ligand_atoms_xyz, params_dict)
        counter = 0
        non_zero = 0
        non_clash = 0
        for u, point in enumerate(binding_site_grid):
            for t, lig_atoms_rotated in enumerate(ligand_orientations):
                cf = get_cf(point, lig_atoms_rotated, energy_matrix, target_grid, min_xyz, cell_width, target_atoms_xyz,
                            pred_05, pred_95, ligand_atoms_numbers_types, target_atoms_numbers_types,
                            atoms_radius, ligand_atoms_radius)
                cfs_list[counter][0] = cf
                cfs_list[counter][1] = t
                cfs_list[counter][2] = u
                counter += 1
        for cf in cfs_list:
            if cf[0] != 0 and cf[0] != 1000000.0:
                non_zero += 1
            if cf[0] != 1000000.0:
                non_clash += 1
        total_evals_list.append(len(cfs_list))
        non_zero_list.append(non_zero)
        non_clash_list.append(non_clash)
        if verbose and Time:
            print(f"REMARK time one ligand: {datetime.datetime.now().time()}")
        index = np.argmin(cfs_list, axis=0)
        cfs_list_by_ligand[i] = cfs_list[index[0]][0]
        if params_dict["OUTPUT_PDB"] == "True":
            translated_coords = np.zeros((len(ligand_orientations[int(cfs_list[index[0]][1])]), 3), dtype=np.float32)
            for atom in range(len(ligand_orientations[int(cfs_list[index[0]][1])])):
                translated_coords[atom] = np.add(ligand_orientations[int(cfs_list[index[0]][1])][atom], binding_site_grid[int(cfs_list[index[0]][2])])
            write_test(translated_coords, name_list[i], './ligand_poses/')

    if verbose and Time:
        print("REMARK time end:")
        print(datetime.datetime.now().time())
    for a, ligand in enumerate(ligands_list):
        ligand_list[a] = ligand.replace("./test_ligand/", "")
    ##############################################################################################
    print(f"REMARK | {'Name':^20} | {'CF':^20} | {'Atoms':^5} | {'Evals no clash and cf != 0':^26} | {'Evals no clash':^14} | {'Total evals':^11}")
    for z, ligand in enumerate(ligand_list):
        print(f"RESULT | {name_list[z]:^20} | {cfs_list_by_ligand[z]:^20} | {atm_quantity[z]:^5} | {non_zero_list[z]:^26} | {non_clash_list[z]:^14} | {total_evals_list[z]:^11}")
    ##############################################################################################


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    verbose = True
    category = None
    Time = False
    atm_quantity = []

    path_to_ligands = sys.argv[1]
    ligand_list = sys.argv[4].split(',')
    name_list = sys.argv[6].split(",")
    receptor = sys.argv[2]
    binding_site = sys.argv[3]
    last_line = sys.argv[5]

    if path_to_ligands.find("active") != -1:
        category = "active"
    if path_to_ligands.find("decoy") != -1:
        category = "decoy"

    config_file = "./config.txt"

    #sys.stdout = open(f"./results/{path_to_ligands.split('/')[-2]}/{category}_{ligand_list[0].replace(path_to_ligands, '')}.txt", "w")
    if verbose and Time:
        print("REMARK time start:")
        print(datetime.datetime.now().time())
    print(f"REMARK Receptor path: {receptor}")
    print(f"REMARK Binding site grid path: {binding_site}")
    main(config_file, ligand_list, binding_site, receptor, './MC_5p_norm_P10_M2_2.txt', last_line)
    #sys.stdout.close()
