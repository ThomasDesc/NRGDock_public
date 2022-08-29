import os
import shutil
import sys

import fix_atom_pdb
from analyse_new_results import main as ana


class RESULT:
    def __init__(self, name, cf, atoms, eval1, eval2, act_or_dec):
        self.name = name
        self.cf = float(cf)
        self.atoms = int(atoms)
        self.eval1 = int(eval1)
        self.eval2 = int(eval2)
        self.act_or_dec = act_or_dec


def get_name_list(target_path):
    ligand_list = []
    info_list = []
    filenames = next(os.walk(target_path), (None, None, []))[2]
    for a, filename in enumerate(filenames):
        result_list = []
        with open(target_path + filename) as f:
            lines = f.readlines()
            for line in lines:
                if a == 0 and line.startswith("REMARK"):
                    info_list.append(line)
                elif line.startswith("RESULT"):
                    itemised_list = line.strip().replace(" ", "").split("|")[1:]
                    result_list.append(itemised_list)
            if filename.split("_")[0] == "active":
                for element in result_list:
                    ligand_list.append(RESULT(element[0], element[1], element[2], element[3], element[4], "active"))
            elif filename.split("_")[0] == "decoy":
                for element in result_list:
                    ligand_list.append(RESULT(element[0], element[1], element[2], element[3], element[4],  "decoy"))
            else:
                ligand_list.append(RESULT(element[0], element[1], element[2], element[3], element[4], "ligand"))
    return ligand_list, info_list


def compress(output_path, ligand_list, target, info_list):
    final_path = f"{output_path}{target}.txt"
    ligand_list.sort(key=lambda x: x.cf)
    with open(final_path, "w") as h:
        for item in info_list:
            h.write(item)
        for ligand in ligand_list:
            h.write(f"RESULT | {ligand.name:^20} | {ligand.cf:^20} | {ligand.atoms:^5} | {ligand.eval1:^26} | {ligand.eval2:^14} | {ligand.act_or_dec:^11}\n")


def check_result_folders_existence(path, tgt):
    if not os.path.exists(path + tgt):
        os.makedirs(path + tgt)


def get_output_name(cfg_path):
    params_dict = dict()
    with open(cfg_path) as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('N_ORIENTATIONS'):
            params_dict['n_orientations'] = str(line.split()[-1])
        elif line.startswith('DOT_DIVISION'):
            params_dict['DOT_DIVISION'] = str(line.split()[-1])
        elif line.startswith('KEPT_PDB_NUMBER'):
            params_dict['KEPT_PDB_NUMBER'] = line.split()[-1]
        elif line.startswith('CLEAN'):
            params_dict['CLEAN'] = line.split()[-1]
    output_path = f"./results_processed/{params_dict['n_orientations']}_rotations_{params_dict['DOT_DIVISION']}_grid/"

    if os.path.exists(output_path):
        number = 2
        while os.path.exists(f"{output_path[:-1]}_{str(number)}/"):
            number += 1
        new_path = f"{output_path[:-1]}_{str(number)}/"
        os.makedirs(new_path)
        output_path = new_path
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    return output_path, params_dict


def reset_result_folder(results, params_dict):
    folder_names = next(os.walk(results), (None, None, []))[1]
    if params_dict["CLEAN"] == 'True':
        if os.path.exists(results):
            filenames = next(os.walk(results), (None, None, []))[1]
            for file in filenames:
                shutil.rmtree(os.path.join(results, file))
    # PREPARE OUTPUT FOLDER
    for folder in folder_names:
        check_result_folders_existence(results, folder)
    # END PREP


def get_good_ligands(number, ligand_list):
    good_ligands = []
    for ligand_counter, ligand in enumerate(ligand_list):
        if ligand_counter > int(number) - 1:
            return good_ligands
        else:
            good_ligands.append(ligand.name)


def delete_ligands(number, ligand_list, output_path, dir):
    good_ligands = get_good_ligands(number, ligand_list)
    new_dir = os.path.join(output_path, dir + "_ligand_poses")
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)
    for filename in good_ligands:
        if os.path.exists(f"./ligand_poses/{dir}"):
            init_dir = f"./ligand_poses/{dir}/"
        else:
            init_dir = f"./ligand_poses/"
        try:
            shutil.copyfile(f"{init_dir}{filename}.pdb", f"{new_dir}/{filename}.pdb")
        except:
            print("error copying ligand")
            continue
        if filename.startswith("CHEMBL"):
            fix_atom_pdb.main(f"{init_dir}/", f"../{dir}/actives_final.mol2", filename)
        elif filename.startswith("ZINC"):
            fix_atom_pdb.main(f"{init_dir}/", f"../{dir}/decoys_final.mol2", filename)
    shutil.rmtree("./ligand_poses")
    os.mkdir("./ligand_poses")


def main(target, result_path, config_file):

    output_path, params_dict = get_output_name(config_file)
    print("output path: ", output_path)
    if target:
        dirs = [target]
    else:
        dirs = os.listdir(result_path)
    for dir in dirs:
        ligand_list, info_list = get_name_list(result_path + dir + "/")
        compress(output_path, ligand_list, dir, info_list)
        delete_ligands(params_dict["KEPT_PDB_NUMBER"], ligand_list, output_path, dir)
    reset_result_folder(result_path, params_dict)
    ana(target, output_path, divide_by_atm_count=False)



if __name__ == "__main__":
    try:
        target = str(sys.argv[1])
    except:
        print('No target specified. Analysing all targets')
        target = None

    result_path = "./results/"
    config_file = "config.txt"

    main(target, result_path, config_file)

