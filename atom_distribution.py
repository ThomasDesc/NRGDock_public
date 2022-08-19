import matplotlib.pyplot as plt
import statistics
from scipy.stats import norm

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
        if line.startswith('@<TRIPOS>BOND') or line.startswith("@<TRIPOS>UNITY"):
            coord_start = 0
            atom_number_list.append(n_atoms)
        if coord_start == 1 and line.startswith('@<TRIPOS>ATOM') == False and line[8] != 'H':
            n_atoms += 1
    return atom_number_list


def count_active_and_decoy(path):
    counter = 0
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        if line == '@<TRIPOS>MOLECULE\n':
            counter += 1
    print(counter)


if __name__ == "__main__":
    active_path = "/home/thomas/Desktop/diverse/cp3a4/actives_final.mol2"
    decoy_path = "/home/thomas/Desktop/diverse/cp3a4/decoys_final.mol2"

    active_number_list = return_list_atom_number(active_path)
    decoy_number_list = return_list_atom_number(decoy_path)

    count_active_and_decoy(active_path)
    count_active_and_decoy(decoy_path)


    mean_decoy = statistics.mean(decoy_number_list)
    mean_active = statistics.mean(active_number_list)
    sd_decoy = statistics.stdev(decoy_number_list)
    sd_active = statistics.stdev(active_number_list)
    active_number_list = sorted(active_number_list)
    decoy_number_list = sorted(decoy_number_list)
    #plt.hist(decoy_number_list, bins=125, label="decoy")
    #plt.hist(active_number_list, bins=125, label="active")
    plt.plot(active_number_list, norm.pdf(active_number_list, mean_active, sd_active), label="active")
    plt.plot(decoy_number_list, norm.pdf(decoy_number_list, mean_decoy, sd_decoy), label="decoy")
    plt.ylim(0, 0.1)
    plt.xlabel('Number of atoms')
    plt.ylabel('Frequency')
    plt.title("Original ligands (from DUD-E)")
    plt.legend()
    plt.tight_layout()
    plt.show()