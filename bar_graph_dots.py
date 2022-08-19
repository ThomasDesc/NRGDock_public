import os

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def use_struct(text_file):
    list = []
    dict = {}
    with open(text_file) as a:
        lines = a.readlines()
        for line in lines:
            if line.startswith("Target"):
                target = line.split(" ")[-1].strip()
            if line.startswith("EF"):
                EF = float(line.split(" ")[-1].strip())
                if EF == 0:
                    EF = 0.1
                dict["Name"] = target
                dict["EF"] = EF
                list.append(dict)
                dict = {}
    return list



def make_bargraph():
    result_path = "/home/thomas/Desktop/New_binding_software/results_processed/orientation_x_dot/"
    name_list = []
    divisions = next(os.walk(result_path))[1]
    list = [[0.01] * len(divisions) for _ in range(8)]
    for a, division in enumerate(divisions):
        analysis_path = result_path + division + "/analysis_results.txt"
        dictionnary = use_struct(analysis_path)
        dict_sorted = sorted(dictionnary, key=lambda x: x.get('Name'))
        for target, item in enumerate(dict_sorted):
            list[target][a] = item["EF"]
    for c, item in enumerate(dict_sorted):
        list[c] = [item["Name"]] + list[c]
    for d, name in enumerate(divisions):
        divisions[d] = name.replace("_", " ")
    divisions = ["Name"] + divisions
    df = pd.DataFrame(list, columns=divisions)
    df.to_excel("/home/thomas/Desktop/output.xlsx")
    print(list)
    print(list[0])
    # set width of bar
    #df.plot(x="Name", y=['0.5 A', '0.75 A', '1 A', '1.25 A', '1.5 A baseline'], kind="bar", figsize=(8, 4), width=0.6)
    #df.plot(x="Name", y=['2 orientation', '4 orientation', '6 orientation', '12 orientation', '18 orientation'], kind="bar", figsize=(8, 4), width=0.6)
    df.plot(x="Name", y=['0.75 x 2', "0.75 x 6", '0.75 x 12'], kind="bar", figsize=(8, 4), width=0.6)
    plt.ylim(0, 7.5)
    plt.xlabel("Target")
    plt.ylabel("EF")
    plt.tight_layout()
    plt.legend(loc='upper right')
    #plt.savefig('/home/thomas/Documents/evals_with_clash/' + target + '.png')
    plt.show()


if __name__ == "__main__":
    make_bargraph()

