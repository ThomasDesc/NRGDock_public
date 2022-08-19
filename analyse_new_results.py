import os
import sys

# import matplotlib.pyplot as plt
# import statistics
# from scipy.stats import norm


def make_unprocessed_list(path, divide_by_atm_count):
    unprocessed_list = []
    active_counter = 0
    decoy_counter = 0
    with open(path) as f:
        lines = f.readlines()
        for line in lines:
            if not line.startswith("REMARK"):
                line = line.strip()
                if line.startswith("RESULT"):
                    line = line.replace(" ","").split('|')[1:]
                    name = line[0]
                    cf = line[1]
                    atm = line[2]
                    if divide_by_atm_count:
                        cf = float(cf) / int(atm)
                    else:
                        cf = float(cf)
                    if line[-1] == "active":
                        unprocessed_list.append(["active_" + name, float(cf)])
                        active_counter += 1
                    if line[-1] == "decoy":
                        unprocessed_list.append(["decoy_" + name, float(cf)])
                        decoy_counter += 1
    return unprocessed_list, active_counter, decoy_counter


def make_ligand_list(unprocessed_list):
    ligand_list = []
    for ligand in unprocessed_list:
        dict = {}
        dict['Name'] = ligand[0]
        dict['CF'] = ligand[1]
        ligand_list.append(dict)
    ligand_list_sorted = sorted(ligand_list, key=lambda x: float(x.get('CF')))
    return ligand_list_sorted


def calculate_EF(ligand_list_sorted, target, active_number, decoy_number):
    active_counter = 0
    decoy_counter = 0
    EF_percentage = 0.01
    organised_results = []
    dict = {}
    for x, ligand in enumerate(ligand_list_sorted):
        if ligand["Name"].startswith("decoy"):
            decoy_counter += 1
        if ligand["Name"].startswith("active"):
            active_counter += 1
        if decoy_counter > len(ligand_list_sorted) * EF_percentage:
            break
    proportion = active_counter/(len(ligand_list_sorted)*EF_percentage)
    EF = str(proportion/(active_number/decoy_number))
    print("EF =", EF)
    print(f"EF % is: {EF_percentage*100}%")
    dict["Name"] = target.replace(".txt", "")
    dict["EF"] = EF
    organised_results.append(dict)
    return organised_results


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
    list_of_dict = [item for sublist in list_of_dict for item in sublist]
    organised_results_sorted = sorted(list_of_dict, key=lambda x: x.get('Name'))
    with open(result_save_path + "analysis_results.txt", "w") as f:
        for result in organised_results_sorted:
            f.write("---------------------------------------")
            f.write("\nTarget = " + result["Name"])
            f.write("\nEF = {:0.2f}\n".format(float(result["EF"])))


def main(target, result_path, divide_by_atm_count):
    organised_results = []
    if target:
        dirs = [f"{target}.txt"]
    else:
        dirs = next(os.walk(result_path))[1]
    for dir in dirs:
        base_path = f"{result_path}{dir}"
        print("\n---------------------------------------")
        print("Target = " + dir)
        processed_list, active_counter, decoy_counter = make_unprocessed_list(base_path, divide_by_atm_count)
        ligand_list_sorted = make_ligand_list(processed_list)
        ef = calculate_EF(ligand_list_sorted, dir, active_counter, decoy_counter)
        organised_results.append(ef)
        # graph(ligand_list_sorted)
    save_results(organised_results, result_path)


if __name__ == "__main__":
    try:
        target = str(sys.argv[1])
    except:
        print('No target specified. Analysing all targets')
        target = None
    divide_by_atm = False
    results = "./results_processed/8_orientations_1.5_grid_2/"
    main(target, results, divide_by_atm)
