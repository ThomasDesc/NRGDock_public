from os import path

def read_ligand_to_fix(name, lig_path):
    with open(path.join(lig_path, name)) as f:
        lines = f.readlines()
    return lines


def find_good_atm_names(name, big_path):
    final_list = []
    if len(name.split("_")) > 1:
        name = name.split("_")[0]
    with open(big_path) as f:
        lines = f.readlines()
        for a, line in enumerate(lines):
            if line.startswith(name):
                atom_lines = lines[a:]
                break

        for b, line in enumerate(atom_lines):
            if line.startswith("@<TRIPOS>ATOM"):
                start = b + 1
            if line.startswith('@<TRIPOS>BOND') or line.startswith("@<TRIPOS>UNITY"):
                end = b
                break
        info = atom_lines[start:end]
        for line in info:
            if line[8] != 'H':
                final_list.append(line)
        return final_list


def execute_fix(to_be_fixed, good_list, output_path):
    fixed_list = []
    for a, line in enumerate(to_be_fixed):
        atom = good_list[a][46:52].replace(" ", "").split(".")[0]
        while len(atom) < 4:
            atom += " "
        corrected_line = line[0:13] + atom + line[17:]
        fixed_list.append(corrected_line)
    with open(output_path, "w") as f:
        for line in fixed_list:
            f.write(line)

def main(path_ligands, path_good_files, name):
    print("fix")
    to_be_fixed = read_ligand_to_fix(name + ".pdb", path_ligands)
    good_list = find_good_atm_names(name, path_good_files)
    execute_fix(to_be_fixed, good_list, path.join(path_ligands, name +".pdb"))




if __name__ == "__main__":
    main("C:/Users/thoma/Desktop/8_orientations_1.5_grid/akt1_ligand_poses", r"C:\Users\thoma\Desktop\diverse\akt1\actives_final.mol2", "CHEMBL211018")
