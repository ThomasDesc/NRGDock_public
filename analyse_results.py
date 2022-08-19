import os
import sys
# import matplotlib.pyplot as plt
# import statistics
# from scipy.stats import norm
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


def calculate_EF(ligand_list_sorted, target):
    active_counter = 0
    decoy_counter = 0
    EF_percentage = 0.01
    dict = {}
    for x, ligand in enumerate(ligand_list_sorted):
        if ligand["Name"].startswith("ZINC"):
            decoy_counter += 1
        if ligand["Name"].startswith("CHEMBL"):
            active_counter += 1
        if decoy_counter > decoy_ligand_quantity * EF_percentage:
            break
    proportion = active_counter/(len(ligand_list_sorted)*EF_percentage)
    EF = str(proportion/(active_ligand_quantity/decoy_ligand_quantity))
    print("EF =", EF)
    print("EF % is: {}%".format(EF_percentage*100))
    dict["Name"] = target
    dict["EF"] = EF
    organised_results.append(dict)


'''def graph(ligand_list_sorted):
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
    plt.show()'''


def save_results(list_of_dict, result_save_path):
    organised_results_sorted = sorted(list_of_dict, key=lambda x: x.get('Name'))
    with open(result_save_path + "analysis_results.txt", "w") as f:
        for result in organised_results_sorted:
            f.write("---------------------------------------")
            f.write("\nTarget = " + result["Name"])
            f.write("\nEF = {:0.2f}\n".format(float(result["EF"])))

if __name__ == "__main__":
    try:
        target = str(sys.argv[1])
    except:
        print('No target specified. Analysing all targets')
        target = None
    result_path = "../results_processed/8_orientations_2.5_grid/"
    organised_results = []
    if target:
        dirs = [target]
    else:
        dirs = next(os.walk(result_path))[1]
    devide_by_atm_count = False
    for dir in dirs:
        base_path = result_path + dir + "/"
        active_path = base_path + "active.txt"
        decoy_path = base_path + "decoy.txt"
        print("\n---------------------------------------")
        print("Target = " + dir)
        active_unprocessed_list = make_unprocessed_list(active_path, "CHEMBL_")
        decoy_unprocessed_list = make_unprocessed_list(decoy_path, "ZINC_")

        active_ligand_quantity = len(active_unprocessed_list)
        decoy_ligand_quantity = len(decoy_unprocessed_list)
        merged_unprocessed_list = active_unprocessed_list + decoy_unprocessed_list

        ligand_list_sorted = make_ligand_list(merged_unprocessed_list)
        calculate_EF(ligand_list_sorted, dir)

        #graph(ligand_list_sorted)
    save_results(organised_results, result_path)