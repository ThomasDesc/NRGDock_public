import numpy as np
from numba import njit

# def njit(f):
#     return f

@njit
def get_surface(distance, energy, pred_05, pred_95):
    index = (distance * 100000).astype(np.int32)[0]
    if energy > 0:
        surface = pred_05[index]
    else:
        surface = pred_95[index]
    return surface


@njit
def get_emat_value(atom_a, atom_b, emat):
    energy_value = emat[int(atom_a)][int(atom_b)]
    return energy_value


@njit
def reset_ligand_coords(ligand_atoms_rotated, grid_point):
    for i in range(len(ligand_atoms_rotated)):
        ligand_atoms_rotated[i] -= grid_point

@njit
def get_cf(grid_point, ligand_atoms_rotated, energy_matrix, target_grid, min_xyz, cell_width, target_atoms_xyz, pred_05,
           pred_95, ligand_atoms_numbers_types, target_atoms_numbers_types, atoms_radius, ligand_atoms_radius,
           water_radius=1.4):
    for i in range(len(ligand_atoms_rotated)):
        ligand_atoms_rotated[i] += grid_point
    cf = 0.0
    for a, atom in enumerate(ligand_atoms_rotated):
        grid_index = ((atom - min_xyz) / cell_width).astype(np.int32)
        for i_offset in [-1, 0, 1]:
            for j_offset in [-1, 0, 1]:
                for k_offset in [-1, 0, 1]:
                    i = i_offset + grid_index[0]
                    j = j_offset + grid_index[1]
                    k = k_offset + grid_index[2]
                    if i < len(target_grid) and j < len(target_grid[0]) and k < len(target_grid[0][0]):
                        if i >= 0 and j >= 0 and k >= 0:
                            if target_grid[i][j][k][0] == -1:
                                continue
                            else:
                                for neighbour in target_grid[i][j][k]:
                                    if neighbour == -1:
                                        break
                                    else:
                                        dist = np.linalg.norm(target_atoms_xyz[neighbour] - atom)
                                        if dist <= 2:
                                            cf = 1000000.0
                                            reset_ligand_coords(ligand_atoms_rotated, grid_point)
                                            return cf
                                        energy_value = get_emat_value(ligand_atoms_numbers_types[a][1], target_atoms_numbers_types[neighbour][1], energy_matrix)
                                        normalised_dist = dist / (atoms_radius[neighbour] + ligand_atoms_radius[a] + (2*water_radius))
                                        if normalised_dist > 1:
                                            continue
                                        else:
                                            cf += energy_value * get_surface(normalised_dist, energy_value, pred_05, pred_95)
    reset_ligand_coords(ligand_atoms_rotated, grid_point)
    return cf