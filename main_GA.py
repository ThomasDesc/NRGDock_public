import numpy as np
from numba import njit
import os
from nrgten.atom import Atom
from complementarity_function import get_cf
import math as m
import datetime
import sys
import main as main_py


def njit(f):
    return f


@njit
def rotate_ligand(ligand_atoms_xyz, params):
    ligand_atoms_xyz = main_py.center_coords(ligand_atoms_xyz)
    divisions = params['N_ORIENTATIONS']
    single_rotation = 360/divisions
    rotation_counter = 0
    ligand_orientations = np.zeros((divisions, divisions, divisions, len(ligand_atoms_xyz), 3), dtype=np.float64)
    for a, x in enumerate(range(divisions)):
        for b, y in enumerate(range(divisions)):
            for c, z in enumerate(range(divisions)):
                x = m.radians(x*single_rotation)
                y = m.radians(y*single_rotation)
                z = m.radians(z*single_rotation)
                rotation_matrix = main_py.Rx(x) * main_py.Ry(y) * main_py.Rz(z)
                for i, coord in enumerate(ligand_atoms_xyz):
                    ligand_orientations[a][b][c][i] = np.array(np.dot(rotation_matrix, coord))
                #write_test(ligand_orientations[rotation_counter], "ligand_rotation_" + str(rotation_counter))
                rotation_counter += 1
    return ligand_orientations


def gen_new_gen(best_cfs):
    new_generation = np.zeros((len(best_cfs)*2, 6), dtype=np.float64)
    for a, individual in enumerate(best_cfs):
        index = a*2
        new_generation[index] = individual[1:]
        delta_indices = [int(np.round(x, 0)) for x in np.random.normal(0, np.abs(np.random.normal(1, 0.5)), 3)]
        new_generation[index+1][0:3] = individual[1:4]

        new_generation[index+1][3] = individual[4] + delta_indices[0]
        new_generation[index+1][4] = individual[5] + delta_indices[1]
        new_generation[index+1][5] = individual[6] + delta_indices[2]
        for b, item in enumerate(new_generation[index+1][3:]):
            if item >= params_dict["N_ORIENTATIONS"]:
                new_generation[index+1][3+b] = params_dict["N_ORIENTATIONS"] - 1
            elif item < 0:
                new_generation[index+1][3+b] = 0
    #TODO: problem here new gen issue
    return new_generation

def calc_cf_sort(binding_site_grid, ligand_orientations, energy_matrix, target_grid, min_xyz, cell_width,
                target_atoms_xyz, pred_05, pred_95, ligand_atoms_numbers_types, target_atoms_numbers_types,
                atoms_radius, ligand_atoms_radius):
    n_cf_evals = len(binding_site_grid) * params_dict["N_ORIENTATIONS"] ** 3
    cfs_list = np.zeros((n_cf_evals, 7), dtype=np.float64)
    counter = 0
    for point in binding_site_grid:
        for a, _ in enumerate(ligand_orientations):
            for b, _ in enumerate(ligand_orientations[a]):
                for c, _ in enumerate(ligand_orientations[a][b]):
                    cf = get_cf(point, ligand_orientations[a][b][c], energy_matrix, target_grid, min_xyz, cell_width,
                                target_atoms_xyz, pred_05, pred_95, ligand_atoms_numbers_types,
                                target_atoms_numbers_types, atoms_radius, ligand_atoms_radius)
                    cfs_list[counter][0] = cf
                    cfs_list[counter][1] = point[0]
                    cfs_list[counter][2] = point[1]
                    cfs_list[counter][3] = point[2]
                    cfs_list[counter][4] = a
                    cfs_list[counter][5] = b
                    cfs_list[counter][6] = c
                    counter += 1
    sorted_cfs_list = cfs_list[cfs_list[:, 0].argsort()]
    top_10_index = 30
    best_cfs = sorted_cfs_list[:top_10_index]

    return best_cfs


def ga_second_pass(binding_site_grid, best_cfs, ligand_orientations, energy_matrix, target_grid, min_xyz, cell_width,
                   target_atoms_xyz, pred_05, pred_95, ligand_atoms_numbers_types, target_atoms_numbers_types,
                   atoms_radius, ligand_atoms_radius, PATH):
    new_generation = gen_new_gen(best_cfs)
    n_cf_evals = len(binding_site_grid) * len(new_generation)
    cfs_list = np.zeros((n_cf_evals, 7), dtype=np.float64)
    counter = 0
    print("\n----------New generation----------")
    print("CF from last generation: ", best_cfs[0][0])
    for point in binding_site_grid:
        for individual in new_generation:
            cf = get_cf(point, ligand_orientations[int(individual[3])][int(individual[4])][int(individual[5])], energy_matrix, target_grid, min_xyz, cell_width,
                        target_atoms_xyz, pred_05, pred_95, ligand_atoms_numbers_types,
                        target_atoms_numbers_types, atoms_radius, ligand_atoms_radius)
            cfs_list[counter][0] = cf
            cfs_list[counter][1] = point[0]
            cfs_list[counter][2] = point[1]
            cfs_list[counter][3] = point[2]
            cfs_list[counter][4:] = individual[3:]
            counter += 1
    sorted_cfs_list = cfs_list[cfs_list[:, 0].argsort()]
    top_10_index = 50
    best_cfs = sorted_cfs_list[:top_10_index]
    # TODO: fix
    for i in range(len(ligand_orientations[int(best_cfs[0][4])][int(best_cfs[0][5])][int(best_cfs[0][6])])):
        ligand_orientations[int(best_cfs[0][4])][int(best_cfs[0][5])][int(best_cfs[0][6])][i] += best_cfs[0][1:4]
    main_py.write_test(ligand_orientations[int(best_cfs[0][4])][int(best_cfs[0][5])][int(best_cfs[0][6])], "best", PATH)
    for i in range(len(ligand_orientations[int(best_cfs[0][4])][int(best_cfs[0][5])][int(best_cfs[0][6])])):
        ligand_orientations[int(best_cfs[0][4])][int(best_cfs[0][5])][int(best_cfs[0][6])][i] -= best_cfs[0][1:4]
    print("\nCF from current generation: ", best_cfs[0][0])
    return best_cfs


def generate_new_dots(best_cfs, ga_sphere_radius):
    sphere_dots = np.zeros((len(best_cfs) * 7, 3), dtype=np.float32)
    for x, dot in enumerate(best_cfs):
        x *= 7
        array = np.random.normal(loc=0, scale=ga_sphere_radius, size=[6, 3])
        for a in range(0, 6):
            sphere_dots[x + a] = np.add(dot[1:4], array[a])
        sphere_dots[x + 6][0] = dot[1]
        sphere_dots[x + 6][1] = dot[2]
        sphere_dots[x + 6][2] = dot[3]
    main_py.write_test(sphere_dots, "sphere_dots", PATH)
    return sphere_dots


@njit
def main(config_file, ligands_list, binding_site, target, energy_matrix_file, last_line):
    params_dict = main_py.get_params_dict(config_file)
    print("REMARK software: main_GA")
    print("REMARK  orientations: ", params_dict["N_ORIENTATIONS"])
    print("REMARK  dot separation: ", params_dict["DOT_DIVISION"])
    print("REMARK  GA generations: ", params_dict["GA_GENERATIONS"])
    print("REMARK  GA sphere radius: ", params_dict["GA_SPHERE_RADIUS"])
    rad_dict = main_py.load_rad_list(PATH + "FastAID_Py/radius_list.txt")
    pred_05 = main_py.import_pred_list("pred_05.txt", params_dict, PATH)
    pred_95 = main_py.import_pred_list("pred_05.txt", params_dict, PATH)
    target_atoms_xyz, target_atoms_numbers_types, atoms_radius = main_py.load_atoms_mol2(target, None, None, rad_dict)
    target_grid, min_xyz, cell_width = main_py.build_3d_cube_grid(params_dict, target_atoms_xyz, atoms_radius)
    OG_binding_site_grid = main_py.load_binding_site_grid(params_dict, binding_site, target_grid, PATH)
    binding_site_grid = main_py.clean_bindig_site_grid(params_dict, target_grid, OG_binding_site_grid,
                                                       min_xyz, cell_width, target_atoms_xyz, verbose, PATH)
    energy_matrix = main_py.load_energy_matrix(energy_matrix_file)
    cfs_list_by_ligand = np.zeros((len(ligands_list)), dtype=np.float32)
    total_evals_list = []
    ga_sphere_radius = params_dict["GA_SPHERE_RADIUS"]
    if verbose and Time:
        print("REMARK ligand start:")
        print(datetime.datetime.now().time())

    for i, ligand in enumerate(ligands_list):
        temp_evals = 0
        ligand_atoms_xyz, ligand_atoms_numbers_types, ligand_atoms_radius = main_py.load_atoms_mol2(path_to_ligands,
                                                                                                    ligand, last_line,
                                                                                                    rad_dict)
        atm_quantity.append(len(ligand_atoms_xyz))
        ligand_orientations = rotate_ligand(ligand_atoms_xyz, params_dict)
        #First GA pass
        best_cfs = calc_cf_sort(binding_site_grid, ligand_orientations, energy_matrix, target_grid, min_xyz,
                                cell_width, target_atoms_xyz, pred_05, pred_95, ligand_atoms_numbers_types,
                                target_atoms_numbers_types, atoms_radius, ligand_atoms_radius)
        temp_evals += len(binding_site_grid) * params_dict["N_ORIENTATIONS"] ** 3
        #main_py.write_test(binding_site_grid, "starting_grid", PATH)
        for g in range(0, params_dict["GA_GENERATIONS"]):
            sphere_dots = generate_new_dots(best_cfs, ga_sphere_radius)
            best_cfs = ga_second_pass(sphere_dots, best_cfs, ligand_orientations, energy_matrix, target_grid, min_xyz,
                                      cell_width, target_atoms_xyz, pred_05, pred_95, ligand_atoms_numbers_types,
                                      target_atoms_numbers_types, atoms_radius, ligand_atoms_radius, PATH)
            temp_evals += len(best_cfs)*2 * len(sphere_dots)
        #end GA
        total_evals_list.append(temp_evals)
        if verbose and Time:
            print("REMARK time one ligand:")
            print(datetime.datetime.now().time())
        cfs_list_by_ligand[i] = best_cfs[0][0]
    if verbose and Time:
        print("REMARK time end:")
        print(datetime.datetime.now().time())
    for a, ligand in enumerate(ligands_list):
        ligand_list[a] = ligand.replace(PATH + "FastAID_Py/test_ligand/", "")
    ##############################################################################################
    print(f"REMARK | {'Name':^20} | {'CF':^20} | {'Atoms':^5} | | | {'Total evals':^11}")
    for z, ligand in enumerate(ligand_list):
        print(f"RESULT | {name_list[z]:^20} | {cfs_list_by_ligand[z]:^20} | {atm_quantity[z]:^5} | | | {total_evals_list[z]:^11}")
    ##############################################################################################


if __name__ == "__main__":

    verbose = True
    PATH = None
    config_file = None
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

    config_file = "/home/thomasd/projects/rrg-najmanov/thomasd/New_binding_software/FastAID_Py/config.txt"
    if not os.path.exists(config_file):
        config_file = "/home/thomasd/projects/def-najmanov/thomasd/New_binding_software/FastAID_Py/config.txt"
        if not os.path.exists(config_file):
            config_file = "config.txt"

    device = main_py.get_device(config_file)
    params_dict = main_py.get_params_dict(config_file)

    if device == "windows":
        PATH = params_dict['WINDOWS_PATH']
    if device == "linux":
        PATH = params_dict['LINUX_PATH']
    if device == "narval":
        PATH = params_dict['NARVAL_PATH']
    if device == "beluga":
        PATH = params_dict['BELUGA_PATH']

    sys.stdout = open(PATH + "results/" + path_to_ligands.split("/")[-2] + "/" + category + "_" + ligand_list[0].replace(path_to_ligands, "") + ".txt", "w")
    if verbose and Time:
        print("REMARK time start:")
        print(datetime.datetime.now().time())
    print(f"REMARK Receptor path: {receptor}")
    print(f"REMARK Binding site grid path: {binding_site}")
    main(PATH + "FastAID_Py/config.txt", ligand_list, binding_site, receptor,
         PATH + 'FastAID_Py/MC_5p_norm_P10_M2_2.txt', last_line)
    sys.stdout.close()