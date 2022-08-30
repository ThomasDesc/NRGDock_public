import shutil
import sys
import os

def find_file(end, path):
    filenames = next(os.walk(path), (None, None, []))[2]
    for filename in filenames:
        if filename.endswith(end):
            return filename


def count_molecules(path):
    line_list = []
    ligand_counter = 0
    name_list = []
    same_name_counter = 0

    with open(path) as f:
        lines = f.readlines()
        for t, line in enumerate(lines):
            if line.find("@<TRIPOS>MOLECULE") != -1:
                ligand_counter += 1
                line_list.append(t)
                name = lines[t+1].strip()
                if len(name_list) != 0:
                    previous_name = name_list[-1].split("_")[0]
                    if previous_name != name:
                        name_list.append(name)
                        same_name_counter = 0
                    else:
                        same_name_counter += 1
                        name += "_" + str(same_name_counter)
                        name_list.append(name)
                elif len(name_list) == 0:
                    name_list.append(name)

    return ligand_counter, line_list, name_list


def divisible(ligand_counter, line_list, job_limit, name_list):
    number = 0
    divider = None
    final_name_list = []
    final_list = []

    for a in range(2, ligand_counter):
        if ligand_counter / a < job_limit:
            divider = a
            break
    while number < ligand_counter:
        final_list.append(line_list[number:number + divider + 1])
        final_name_list.append(name_list[number:number + divider])
        number += divider
    return final_list, final_name_list


def build_string_list(final_list, ligand_path, name_list, ga, job_number, receptor, binding_site):
    command_string_list = []
    global next_job_counter
    for a, lig_group in enumerate(final_list):
        string = f"python3 {os.path.join(software_path, 'main_GA.py')} {ligand_path} {receptor} {binding_site} "
        if ga == 'False':
            string = string.replace("main_GA.py", "main.py")
        for c, ligand in enumerate(lig_group[:-1]):
            if c != len(lig_group) - 2:
                string += str(ligand) + ","
            else:
                string += str(ligand)
        string += " " + str(lig_group[-1])  # add line where last ligand .mol ends
        # adding names
        string += " "
        for d, name in enumerate(name_list[a]):
            if d != len(name_list[a]) - 1:
                string += str(name) + ","
            else:
                string += str(name)
        if next_job_counter < job_number:
            string += "\nsbatch job_" + str(next_job_counter) + ".sh"
            next_job_counter += 1
        command_string_list.append(string)
    return command_string_list


def clean_job_folder():
    _, _, files = next(os.walk("./jobs/"), (None, None, []))
    for file in files:
        os.remove("./jobs/" + file)


def build_sbatch_list(command_string_list):
    sbatch_list = []
    d = None
    global job_counter
    for d, command in enumerate(command_string_list):
        with open("./job_template.sh") as f:
            lines = f.readlines()
            with open("./jobs/job_" + str(job_counter + d) + ".sh", "w") as g:
                lines.append(command)
                for line in lines:
                    g.write(line)
                sbatch_list.append("sbatch job_" + str(job_counter + d) + ".sh")
    job_counter += d + 1
    return sbatch_list


def check_output_path_existence(path, target):
    if os.path.exists(path + target):
        shutil.rmtree(path + target)
        os.makedirs(path + target)
    else:
        os.makedirs(path + target)


def change_account(account, software_path):
    with open(os.path.join(software_path, "job_template.sh")) as f:
        lines = f.readlines()
        for a, line in enumerate(lines):
            if line.startswith("#SBATCH --account"):
                lines[a] = f"#SBATCH --account {account}\n"
            if line.startswith("source"):
                path = os.path.dirname(software_path)
                lines[a] = f"source {os.path.join(path,'ENV', 'bin', 'activate')}\n"

    with open(os.path.join(software_path, "job_template.sh"), "w") as g:
        for line in lines:
            g.write(line)


if __name__ == "__main__":

    config_file = "./config.txt"
    software_path = os.path.dirname(os.path.realpath(__file__))
    receptor_path = sys.argv[1]
    if not receptor_path.endswith("/"):
        receptor_path += "/"
    account = sys.argv[2]
    conformer = sys.argv[3]
    GA = sys.argv[4]
    try:
        enriching = sys.argv[5]
    except:
        print("missing argument, True or False, do you have 2 ligand files (active and decoy)?")
        exit()
    next_job_counter = 1000
    if account == "True" or account == "False" or account == "GROUP":
        exit("Incorrect argument. Please input your account name.")
    change_account(account, software_path)

    if conformer == "True":
        active_ligand_path = os.path.join(receptor_path, "actives_final_conf.mol2")
        decoy_ligand_path = os.path.join(receptor_path, "decoys_final_conf.mol2")
    else:
        active_ligand_path = os.path.join(receptor_path, "actives_final.mol2")
        decoy_ligand_path = os.path.join(receptor_path, "decoys_final.mol2")

    # PREPARE OUTPUT FOLDER
    output_path = "./results/"
    ligand_pose_path = "./ligand_poses/"
    check_output_path_existence(output_path, receptor_path.split("/")[-2])
    check_output_path_existence(ligand_pose_path, receptor_path.split("/")[-2])
    # END PREP

    receptor = os.path.join(receptor_path, find_file("receptor.mol2", receptor_path))
    binding_site = os.path.join(receptor_path, "get_cleft", find_file("_sph_1.pdb", f"{receptor_path}/get_cleft/"))

    if enriching == "True":
        active_counter, active_line_list, active_name = count_molecules(active_ligand_path)
        decoy_counter, decoy_line_list, decoy_name = count_molecules(decoy_ligand_path)
        active_job_limit = int(2000 / ((active_counter + decoy_counter) / active_counter))
        decoy_job_limit = int(2000 / ((active_counter + decoy_counter) / decoy_counter)) - 1
        active_final_list, active_final_name_list = divisible(active_counter, active_line_list,
                                                              active_job_limit, active_name)
        decoy_final_list, decoy_final_name_list = divisible(decoy_counter, decoy_line_list,
                                                            decoy_job_limit, decoy_name)
        job_number = len(active_final_list) + len(decoy_final_list)
        active_string_list = build_string_list(active_final_list, active_ligand_path,
                                               active_final_name_list, GA, job_number, receptor, binding_site)
        decoy_string_list = build_string_list(decoy_final_list, decoy_ligand_path,
                                              decoy_final_name_list, GA, job_number, receptor, binding_site)

        clean_job_folder()

        job_counter = 0
        active_sbatch_list = build_sbatch_list(active_string_list)
        decoy_sbatch_list = build_sbatch_list(decoy_string_list)

        total_sbatch_list = active_sbatch_list + decoy_sbatch_list

        with open("./jobs/start_jobs.sh", "w") as h:
            for f, element in enumerate(total_sbatch_list):
                if f < 1000:
                    h.write(element + "\n")
    if enriching == "False":
        active_ligand_path = os.path.join(receptor_path, "ligands.mol2")
        active_counter, active_line_list, active_name = count_molecules(active_ligand_path)
        active_job_limit = 2000
        active_final_list, active_final_name_list = divisible(active_counter, active_line_list,
                                                              active_job_limit, active_name)
        job_number = len(active_final_list)
        active_string_list = build_string_list(active_final_list, active_ligand_path,
                                               active_final_name_list, GA, job_number, receptor, binding_site)

        clean_job_folder()

        job_counter = 0
        active_sbatch_list = build_sbatch_list(active_string_list)

        total_sbatch_list = active_sbatch_list

        with open("./jobs/start_jobs.sh", "w") as h:
            for f, element in enumerate(total_sbatch_list):
                if f < 1000:
                    h.write(element + "\n")
