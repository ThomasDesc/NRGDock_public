import os

import matplotlib.pyplot as plt
import statistics
from scipy.stats import norm


def make_unprocessed_list(path, prefix):
    unprocessed_list = []
    with open(path) as f:
        lines = f.readlines()
        for line in lines:
            if not line.startswith("REMARK"):
                line = line.strip()
                if line.startswith("RESULT"):
                    line = line.replace("RESULT ", "")
                line = line.split(" ")
                if devide_by_atm_count:
                    unprocessed_list.append([prefix + line[0], float(line[1]) / float(line[2])])
                else:
                    unprocessed_list.append([prefix + line[0], float(line[1])])
    return unprocessed_list


def make_ligand_list(unprocessed_list):
    ligand_list = []
    for ligand in unprocessed_list:
        dict = {}
        dict['Name'] = ligand[0]
        dict['CF'] = ligand[1]
        ligand_list.append(dict)
    ligand_list_sorted = sorted(ligand_list, key=lambda x: float(x.get('CF')))
    return ligand_list_sorted


def calculate_EF(ligand_list_sorted):
    active_counter = 0
    decoy_counter = 0
    EF_percentage = 0.01
    dict = {}
    ligand_number_percentage = len(ligand_list_sorted) * EF_percentage
    for x, ligand in enumerate(ligand_list_sorted):
        if ligand["Name"].startswith("ZINC"):
            decoy_counter += 1
        if ligand["Name"].startswith("CHEMBL"):
            active_counter += 1
        if x > ligand_number_percentage:
            break
    proportion = active_counter/(ligand_number_percentage)
    EF = str(proportion/(active_ligand_quantity/decoy_ligand_quantity))
    print("EF =", EF)
    print("EF % is: {}%".format(EF_percentage*100))
    dict["Name"] = "akt1"
    dict["EF"] = EF
    organised_results.append(dict)


def graph(ligand_list_sorted):
    cf_active_list_plot = []
    cf_decoy_list_plot = []
    for ligand in ligand_list_sorted:
        if ligand["Name"].startswith("C"):
            cf_active_list_plot.append(ligand["CF"])
        if ligand["Name"].startswith("Z"):
            cf_decoy_list_plot.append(ligand["CF"])


    mean_decoy = statistics.mean(cf_decoy_list_plot)
    mean_active = statistics.mean(cf_active_list_plot)
    sd_decoy = statistics.stdev(cf_decoy_list_plot)
    sd_active = statistics.stdev(cf_active_list_plot)

    plt.plot(cf_decoy_list_plot, norm.pdf(cf_decoy_list_plot, mean_decoy, sd_decoy), label="decoy")
    plt.plot(cf_active_list_plot, norm.pdf(cf_active_list_plot, mean_active, sd_active), label="decoy")

    plt.xlabel('CF')
    plt.ylabel('Frequency')
    plt.legend(['decoy', 'active'])
    plt.tight_layout()
    plt.show()


def save_results(list_of_dict, result_save_path):
    organised_results_sorted = sorted(list_of_dict, key=lambda x: x.get('Name'))
    with open(result_save_path + "analysis_results.txt", "w") as f:
        for result in organised_results_sorted:
            f.write("---------------------------------------")
            f.write("\nTarget = " + result["Name"])
            f.write("\nEF = {:0.2f}\n".format(float(result["EF"])))

if __name__ == "__main__":
    result_path = "/home/thomas/Desktop/New_binding_software/akt1_test/"
    organised_results = []
    devide_by_atm_count = False
    dirs = next(os.walk(result_path))[1]
    active_path = "/home/thomas/Desktop/New_binding_software/akt1_test/akt1/active.txt"
    active_unprocessed_list = make_unprocessed_list(active_path, "CHEMBL_")
    active_ligand_quantity = len(active_unprocessed_list)
    final_decoy_list = []
    for dir in dirs:
        base_path = result_path + dir + "/"
        decoy_path = base_path + "decoy.txt"
        print("\n---------------------------------------")
        print("Target = " + dir)
        decoy_unprocessed_list = make_unprocessed_list(decoy_path, "ZINC_" + dir + "_")
        final_decoy_list += decoy_unprocessed_list

    decoy_ligand_quantity = len(final_decoy_list)
    merged_unprocessed_list = active_unprocessed_list + final_decoy_list

    ligand_list_sorted = make_ligand_list(merged_unprocessed_list)
    calculate_EF(ligand_list_sorted)

    #graph(ligand_list_sorted)
    save_results(organised_results, result_path)
    None