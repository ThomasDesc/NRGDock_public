import os
import time
from os import walk
import gzip
import shutil
command = []
written_ligands = []


def empty_actives_decoys(folders):
    for folder in folders:
        if os.path.exists(folder + '/actives_split'):
            shutil.rmtree(folder + '/actives_split/')
        if os.path.exists(folder + '/decoys_split'):
            shutil.rmtree(folder + '/decoys_split/')


def writedf(finallist, ligand_name, folder, type):
    counter = 0
    while (ligand_name in written_ligands) == True:
        counter += 1
        ligand_name += '_' + str(counter)
    new_file_path = folder + '/' + type + '/' + ligand_name + '.mol2'
    written_ligands.append(ligand_name)
    with open(new_file_path, "w") as f:
        f.write("".join(finallist))


def untar(folders):
    for folder in folders:
        _, _, files = next(walk(folder), (None, None, []))
        for file in files:
            if file.endswith('.mol2.gz') == True:
                with gzip.open(folder + '/' + file, 'rb') as f_in:
                    if '.gz' in file:
                        extractedfile = file.replace('.gz', '')
                    with open(folder + '/' + extractedfile, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(folder + '/' + file)


def processligand(folders):
    for folder in folders:
        print(folder)
        subfolders = [f.path for f in os.scandir(folder) if f.is_dir()]
        for subfolder in subfolders:
            if subfolder.endswith('actives_split') == True or subfolder.endswith('decoys_split') == True :
                _, _, filenames = next(walk(subfolder), (None, None, []))
                for filename in filenames:
                    if filename.endswith('mol2'):
                        command = "/home/thomas/Desktop/Process_Ligand-master/Process_Ligand -f " + subfolder + '/' + filename +  ' --res_number 9999 --atom_index 90000'
                        os.system(command)


def process_ligand_receptor(folders):
    for folder in folders:
        _, _, filenames = next(walk(folder), (None, None, []))
        for filename in filenames:
            if filename.startswith('receptor_') == True and filename.endswith('.pdb') == True and filename.endswith('inp.pdb') == False:
                run_process_ligand = "/home/thomas/Desktop/Process_Ligand-master/Process_Ligand -f " + folder + '/' + filename + ' -target --res_number 9999 --atom_index 90000'
                os.system(run_process_ligand)


def getcleft(folders):
    for folder in folders:
        if not os.path.exists(folder + '/getcleft'):
            os.makedirs(folder + '/getcleft')
        _, _, filenames = next(walk(folder), (None, None, []))
        for filename in filenames:
            if filename.startswith('receptor_') == True and filename.endswith('.inp.pdb') == False and filename.endswith('.pdb') == True:
                getcleft_command = '/home/thomas/Desktop/Get_Cleft-master/Get_Cleft -p ' + folder + '/' + filename + ' -o ' + folder + '/getcleft/ -s -t 3'
                os.system(getcleft_command)

if __name__ == "__main__":

    folders = [f.path for f in os.scandir('C:/Users/thoma/Desktop/diverse') if f.is_dir()]
    untar(folders)
    empty_actives_decoys(folders)
    for folder in folders:
        written_ligands = []
        print(folder)
        if not os.path.exists(folder +'/actives_split'):
            os.makedirs(folder +'/actives_split')
        if not os.path.exists(folder +'/decoys_split'):
            os.makedirs(folder +'/decoys_split')
        _, _, filenames = next(walk(folder), (None, None, []))

        '''for filename in filenames:
            if filename.find('actives_final.mol2') != -1:
                print("active")
                split_mol2(folder + '/' + filename, folder, 'actives_split')
            if filename.find('decoys_final.mol2') != -1:
                print("decoy")
                split_mol2(folder + '/' + filename, folder, 'decoys_split')'''

    #processligand(folders)
    #process_ligand_receptor(folders)
    #getcleft(folders)
