import time

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import statistics


def average_dots(basepath):
    folders = next(os.walk(base_path))[1]
    for folder in folders:
        next_path = base_path + folder + "/"
        targets = next(os.walk(next_path))[1]
        for target in targets:
            analysis_path = next_path + target + "/decoy.txt"
            list = []
            with open(analysis_path) as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith("REMARK Total"):
                        list.append(int(line.split("  ")[-1].strip()))
        print(folder)
        print(str(statistics.mean(list)) + "\n")


def max_ef(diverse_folder):
    EF_percentage = 0.01
    folders = next(os.walk(diverse_folder))[1]
    for folder in folders:
        active_ligand_count = 0
        decoy_ligand_count = 0
        filenames = next(os.walk(os.path.join(diverse_folder, folder)), (None, None, []))[2]
        for filename in filenames:
            if filename.endswith(".mol2") and filename.startswith("actives"):
                with open(os.path.join(diverse_folder, folder, filename)) as f:
                    lines = f.readlines()
                    for line in lines:
                        if line.find("@<TRIPOS>MOLECULE") != -1:
                            active_ligand_count += 1
            if filename.endswith(".mol2") and filename.startswith("decoys"):
                with open(os.path.join(diverse_folder, folder, filename)) as f:
                    lines = f.readlines()
                    for line in lines:
                        if line.find("@<TRIPOS>MOLECULE") != -1:
                            decoy_ligand_count += 1

        #proportion = active_ligand_count / (active_ligand_count + decoy_ligand_count) #* EF_percentage
        EF = 1 / (active_ligand_count / decoy_ligand_count)

        print("Max EF for target {} is {}".format(folder, EF))

if __name__ == "__main__":
    diverse_folder = "/home/thomas/Desktop/diverse"
    base_path = "/home/thomas/Desktop/New_binding_software/results_processed/grid_dot_variation/"
    #average_dots(base_path)
    max_ef(diverse_folder)