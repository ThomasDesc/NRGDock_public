import time

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

def return_list_atom_number(path):
    atom_number_list = []
    with open(path) as f:
        lines = f.readlines()
    n_atoms = 0
    coord_start = 0
    for line in lines:
        if line.startswith('@<TRIPOS>ATOM'):
            coord_start = 1
            n_atoms = 0
            continue
        if line.startswith('@') and line.startswith('@<TRIPOS>ATOM') == False and coord_start == 1:
            coord_start = 0
            atom_number_list.append(n_atoms)
        if coord_start == 1 and line.startswith('@<TRIPOS>ATOM') == False and line[8] != 'H':
            n_atoms += 1
    return atom_number_list

def make_unprocessed_list(path):
    unprocessed_list = []
    prefix = None
    if path.find("active"):
        prefix = "CHEMBL_"
    else:
        prefix = "ZINC_"
    with open(path) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("REMARK"):
                continue
            else:
                line = line.strip()
                line = line.split(" ")
                unprocessed_list.append([prefix + line[1], float(line[2]), int(line[3]), int(line[4]), int(line[5]), int(line[6])])
    return unprocessed_list

def make_dict_sort_list(unprocessed_list):

    ligand_list = []
    for ligand in unprocessed_list:
        dict = {}
        dict['Name'] = ligand[0]
        dict['CF'] = ligand[1]
        dict["ATOMCOUNT"] = ligand[2]
        dict["NONZERONONCLASH"] = ligand[3]
        dict["NONZERO"] = ligand[4]
        dict["TOTALEVAL"] = ligand[5]
        dict["line_number"] = ligand[0].split("_")[1]
        dict["CLASHES"] = ligand[5] - ligand[4]
        if ligand[5] < ligand[4] :
            print(ligand)
            time.sleep(1)
        ligand_list.append(dict)
    #ligand_list_sorted = sorted(ligand_list, key=lambda x: float(x.get('line_number')))
    #return ligand_list_sorted
    return ligand_list


def count_active_and_decoy(path):
    counter = 0
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        if line == '@<TRIPOS>MOLECULE\n':
            counter += 1
    print(counter)

if __name__ == "__main__":
    #active_path = "/home/thomas/Desktop/diverse/cxcr4/actives_final_conf.mol2"
    #decoy_path = "/home/thomas/Desktop/diverse/cxcr4/decoys_final_conf.mol2"
    base_path = "/home/thomas/Desktop/New_binding_software/results_processed/baseline/"
    targets = next(os.walk(base_path))[1]
    for target in targets:
        active_result_path = "/home/thomas/Desktop/New_binding_software/results_processed/baseline/" + target + "/active.txt"
        decoy_result_path = "/home/thomas/Desktop/New_binding_software/results_processed/baseline/" + target + "/decoy.txt"

        #count_active_and_decoy(active_path)
        #count_active_and_decoy(decoy_path)

        active_number_list = make_unprocessed_list(active_result_path)
        decoy_number_list = make_unprocessed_list(decoy_result_path)

        active_number_list = make_dict_sort_list(active_number_list)
        decoy_number_list = make_dict_sort_list(decoy_number_list)

        total_list = active_number_list + decoy_number_list
        for element in active_number_list:
            element["CF"] /= element["ATOMCOUNT"]
        for element in decoy_number_list:
            element["CF"] /= element["ATOMCOUNT"]
        df = pd.DataFrame(data=total_list)
        ax = sns.violinplot(x="ATOMCOUNT", y="CLASHES", data=df)
        plt.xlabel("Number of atoms in ligand")
        #plt.ylabel("Number of evals with no clash and != 0")
        plt.ylabel("Number of evals with a clash")
        plt.title("Evals with clashes as a function of the # of atoms for " + target)
        plt.tight_layout()
        plt.savefig('/home/thomas/Documents/evals_with_clash/' + target + '.png')
        plt.show()


