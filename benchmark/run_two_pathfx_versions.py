# output_run_two_pathfx_versions.py
# written 12-27-21 AM

from functools import partial
import pandas as pd
import sys
import os
import argparse
from multiprocessing import Pool
import pickle
from random import sample
import shutil
from pathlib import Path
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib import backend_bases
import matplotlib.pyplot as plt


def read_pickle_file(file_path):
    if not os.path.exists(file_path):
        sys.exit("Can't locate input file %s" % file_path)
    return pd.read_pickle(file_path)

def compare_n_drugs(drug_list=None):
    
    # check if a list of drugs was passed as an argument or not
    if drug_list is None or (len(drug_list) == 0):
            print("Drug Input Error: No Drugs given. Please enter valid drug names and DrugBank IDs accepted by both PathFX Versions.\n")
            sys.exit(1)
    num_of_drugs = len(drug_list)
    
    # e.g. paths = ['/Users/aryahaghighi/Documents/PathFX/WorkingBenchmarkPathFX/PathFX_v1/scripts', '/Users/aryahaghighi/Documents/PathFX/WorkingBenchmarkPathFX/PathFx_v2/scripts']
    paths = read_paths()

    # Path to database of all drugs for v1
    # e.g. drug_db_path_v1 = "/Users/aryahaghighi/Documents/PathFX/WorkingBenchmarkPathFX/PathFX_v1/rscs/drug_intome_targets.pkl"
    drug_db_path_v1 = paths[0][0:len(paths) - 10] + "/rscs/drug_intome_targets.pkl"

    # Path to database of all drugs for v2
    # e.g. drug_db_path_v1 = "/Users/aryahaghighi/Documents/PathFX/WorkingBenchmarkPathFX/PathFX_v2/rscs/pfxDB050620_dint.pkl"
    drug_db_path_v2 = paths[1][0:len(paths) - 10] + "/rscs/pfxDB050620_dint.pkl"

    # opening drug-target dictionaries for both PathFX Versions
    drug_db_dict_v1 = read_pickle_file(drug_db_path_v1)
    drug_db_dict_v2 = read_pickle_file(drug_db_path_v2)    
    drug_db_dict_v1_keys = list(drug_db_dict_v1.keys())
    drug_db_dict_v2_keys = list(drug_db_dict_v2.keys())

    # Drug Intersection from both PathFX Versions
    all_unique_drugs_both_pathfx_versions_lst = list(set(drug_db_dict_v1_keys) & set(drug_db_dict_v2_keys))

    # Drugbank IDs from both PathFX Versions
    all_unique_drungbank_ids_both_pathfx_versions_lst = [drug for drug in all_unique_drugs_both_pathfx_versions_lst if drug[0:2] == "DB"]

    print("\nNumber of unique Drugs (names and DrugBank IDs) in both versions of the PathFX interaction databases (NOT LOWERCASED CORRECTED): " + str(len(all_unique_drugs_both_pathfx_versions_lst)))
    print("Number of unique Drugs (DrugBank ID only) in both versions of the PathFX interaction databases (NOT LOWERCASED CORRECTED): " + str(len(all_unique_drungbank_ids_both_pathfx_versions_lst)))
    
    
    """
    
    WARNING!!!!
    
    We used NON-lowercased dictionaries that had less intersections than lowercased dictionaries for both PathFX versions because the original pathFX versions use th non-lowercased dictionaries
    
    """
    
    
    
    # Validate passed in drug list
    for drug in drug_list:
        if drug.lower() not in map(str.lower, all_unique_drugs_both_pathfx_versions_lst):
            print("Drug Input Error: '" + drug + "' is not a valid drug. Please enter valid drug names and DrugBank IDs accepted by both PathFX Versions. \n")
            sys.exit(1)
    
    # Save original working dir before changing it
    ogwd = os.getcwd()

    print("\nOutput from Both PathFX Versions:")
    print("===============================================================================================================================================\n")

    # Run analyses on the randomly selected drugs
    with Pool(20) as pool:
        analysis_v1 = partial(run_analysis_v1, paths)
        analysis_v2 = partial(run_analysis_v2, paths)
        pool.imap_unordered(analysis_v1, drug_list)
        pool.imap_unordered(analysis_v2, drug_list)
        pool.close()
        pool.join()

    print("===============================================================================================================================================\n")

    # Set working dir back to original dir
    os.chdir(ogwd)

    # outputs data to a spreadsheet
    export_data(paths, drug_list)
    
    return None



"""
Checking if the configuration file exists
Reading in pathway to script files for respective PathFX Version
------------------------------------------------------------------------------------
/Users/aryahaghighi/Documents/PathFX/WorkingBenchmarkPathFX/PathFX_v1/scripts
/Users/aryahaghighi/Documents/PathFX/WorkingBenchmarkPathFX/PathFx_v2/scripts
"""
def read_paths():
    locs = ['', '']
    temp_locs = ['', '']
    try:
        with open('pathfx_locs.txt', 'r+') as pathfx_locs:
            locs = pathfx_locs.readlines()
            # If the file is created, but empty
            if len(locs) == 0:
                locs = ['', '']
            # There is only one line in the config file
            if len(locs) == 1:
                locs = [locs[0] + "\n", ""]
            # Throw away all other data besides the first two lines
            locs = [locs[0], locs[1]]
            temp_locs = locs.copy()
            # Check to make sure dirs contain the script file
            for i in range(0, 2):
                path_incorrect = True
                while path_incorrect:
                    try:
                        # check if the script 'phenotype_enrichment_pathway.py' is in the folder of the given pathway
                        if 'phenotype_enrichment_pathway.py' in os.listdir(locs[i].strip()):
                            path_incorrect = False
                        # Input is a path but it doesn't contain the file
                        else:
                            print("The path doesn't have any version of PathFX, try again for v" + i + ": \n")
                            input_loc = input("Enter the path to PathFX version " + i + ": \n")
                            print("\n")
                            input_loc = input_loc.strip()
                            input_loc = input_loc + "/scripts"
                            locs[i] = input_loc + "\n"
                    # Input isn't a path
                    except IOError:
                        input_loc = input("The config doesn't have a valid path for PathFX v" + i + ", enter a path: \n")
                        print("\n")
                        input_loc = input_loc.strip()
                        input_loc = input_loc + "/scripts"
                        locs[i] = input_loc + "\n"
    except IOError:
        # File doesn't exist, so create and then run this function again
        f = open('pathfx_locs.txt', 'w')
        f.close()
        read_paths()
    if locs != temp_locs:
        # Delete current contents and write the new paths
        f = open('pathfx_locs.txt', 'w')
        f.write(locs[0].strip() + "\n" + locs[1])
    # band aid fix for read_paths() returning ["",""] instead of the paths when first creating the config
    if [locs[0].strip(), locs[1].strip()] == ["", ""]:
        return read_paths()
    # Create results folder to hold spreadsheets
    if not Path('results').is_dir():
        os.mkdir('results')
    return [locs[0].strip(), locs[1].strip()]

def run_analysis_v1(paths_in, drug_in):
    # Change working dir so the script can find its resources
    os.chdir(paths_in[0])
    os.system("python3 phenotype_enrichment_pathway.py -d '" + drug_in + "' -a PathFX_analysis")

def run_analysis_v2(paths_in, drug_in):
    # Change working dir so the script can find its resources
    # NOTE: for version 2 of pathfx, file is linked to phenotype_enrichment_pathway_Pfx050120
    os.chdir(paths_in[1])
    os.system("python3 phenotype_enrichment_pathway.py -d '" + drug_in + "' -a PathFX_analysis")


    
    
    
"""
obtaining information from _merged_neighborhood__assoc_table_.txt file for both PathFX versions
and exporting information to an .xlsx for each drug passed by user 
"""
def export_data(input_paths, drug_list):
    
    print("\nExporting '_merged_neighborhood__assoc_table_.txt' files from both PathFX versions to an .xlsx for each drug passed by user\n")
    print("==================================================================")
    for drug in drug_list:
        
        # obtaining file pathways for _merged_neighborhood__assoc_table_ for both PathFX Versions
        paths = ['', '']
        paths[0] = input_paths[0][:-7] + "results/pathfx_analysis/" + drug + "/" + drug + "_merged_neighborhood__assoc_table_.txt"
        paths[1] = input_paths[1][:-7] + "results/pathfx_analysis/" + drug + "/" + drug + "_merged_neighborhood__assoc_table_.txt"

        # Constructing DataFrame for both versions and removing unnecessary columns
        data_frame_v1 = pd.read_csv(paths[0], delimiter="\t", index_col=False)
        data_frame_v2 = pd.read_csv(paths[1], delimiter="\t", index_col=False)
        data_frame_v1.__delitem__('rank')
        data_frame_v2.__delitem__('rank')
        data_frame_v1.__delitem__('Unnamed: 8')
        data_frame_v2.__delitem__('Unnamed: 8')

        # writing Dataframes to the same excel sheet in seperate sheets
        with pd.ExcelWriter("output_run_two_pathfx_versions/" + drug + ".xls") as writer:
            data_frame_v1.to_excel(writer, sheet_name=drug + " all CUI's v1", index=0)
            data_frame_v2.to_excel(writer, sheet_name=drug + " all CUI's v2", index=0)
            
        print("\n\tFinished exporting data to .xlsx file for " + drug)

    print("\nFinished exporting '_merged_neighborhood__assoc_table_.txt' files from both PathFX versions for all drugs passed by user\n")
    print("==================================================================")

    
    
    
    
    
# Setting up argparse
parser = argparse.ArgumentParser(description="Calculate differences between between identified phenotypes and genes to help compare PathFX versions")
command_group_optional_args = parser.add_mutually_exclusive_group()
command_group_optional_args_analysis = parser.add_mutually_exclusive_group()
command_group_optional_args.add_argument('-d', '--drugs', metavar='drug1', type=str, nargs='+', help='List of drugs to compare with PathFX versions. Use quotes if a drug is multiple words. Compares genes by default. Use --phenotypes or --genes-with-phenotype to change this.')

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

if args.drugs:
    compare_n_drugs(drug_list=args.drugs)
print("\nEnd of run_two_pathfx_versions.py script")
print("\n")
















