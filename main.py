import json
import subprocess
import tarfile
from pathlib import Path
import pandas as pd
from src.machine_learning.machine_learning import MachineLearning
from src.conservation_parser import ConservationParser
from src.shape_parser import ShapeParser
from src.rna_folder import RNAFolder

def load_config(path):
    with open("config.json") as config_file:
        config = json.load(config_file)

        settings = config["settings"]
        directories = config["directories"]

        dir_paths = list(directories.values())
        for path in dir_paths:
            Path(path).mkdir(parents = True, exist_ok = True)

        return settings, directories

def determine_mane_version(settings):
    with open("src/ensembl_mane_version_lookup.json") as mane_ver_file:
        ensembl_mane_lookup = json.load(mane_ver_file)
       
        mane_version = "1.0" # default to mane 1.0 if ensembl release specified isn't applicable to try to salvage
        if settings["ensembl_release"] in ensembl_mane_lookup:
            mane_version = ensembl_mane_lookup[settings["ensembl_release"]]
        
        return mane_version

def main():

    print("\nLoading config and setting up...")
    config_path = "config.json"
    settings, directories = load_config(config_path)
    mane_version = determine_mane_version(settings)
    print("Done.\n")


    if settings["use_precompiled_conservation"]:
        print("Readying precompiled conservation data...")
        with tarfile.open("precompiled_conservation_data.tar.gz", "r:gz") as tar:
            tar.extractall(path=".")
        print("Done.\n")

    if settings["use_precompiled_shape"]:
        print("Readying precompiled shape data...")
        with tarfile.open("precompiled_shape_data.tar.gz", "r:gz") as tar:
            tar.extractall(path=".")
        print("Done.\n")


    print("Attempting to download data for annotation...")
    process = subprocess.run(["sh", "./src/download_annotation_data.sh", 
        settings["use_caching"], directories["preload_data"], settings["ensembl_release"], mane_version], shell = False)
    if (process.returncode != 0):
        print("An error occurred while downloading annotation data.")
        return()
    print("Done.\n")


    print("Parsing annotation data...")
    process = subprocess.run(["Rscript", "src/parse_annotation_data.r", config_path], shell = False)
    if (process.returncode != 0):
        print("An error occurred while parsing annotation data.")
        return()
    print("Done.\n")


    print("Extracting 3' UTR and CDS sequences...")
    process = subprocess.run(["Rscript", "src/extract_sequences.r", config_path], shell = False)
    if (process.returncode != 0):
        print("An error occurred while extracting sequences.")
        return()
    print("Done.\n")


    print("Generating conservation scores... This will take a long time if there is no cache.")
    process = subprocess.run(["Rscript", "src/generate_conservation_scores.r", config_path], shell = False)
    if (process.returncode != 0):
        print("An error occurred while generating conservation scores.")
        return()
    print("Done.\n")


    print("Locating binding sites in each experiment...")
    process = subprocess.run(["Rscript", "src/locate_binding_sites.r", config_path], shell = False)
    if (process.returncode != 0):
        print("An error occurred while locating binding sites.")
        return()
    print("Done.\n")


    print("Extracting folding windows for each experiment...")
    process = subprocess.run(["Rscript", "src/extract_windows.r", config_path], shell = False)
    if (process.returncode != 0):
        print("An error occurred while extracting folding windows.")
        return()
    print("Done.\n")


    # print("Folding sequences using ViennaRNA...")
    # RNAFolder(settings, directories)
    # print("Done.\n")


    # print("Extracting features for each experiment...")
    # exit_code = subprocess.call("Rscript src/extract_features.r " + config_path, shell = True)
    # if (exit_code > 0):
    #     print("An error occurred while extracting features.")
    #     return()
    # print("Done.\n")


    # print("Parsing conservation scores for each experiment...")
    # ConservationParser(settings, directories)
    # print("Done.\n")


    # print("Parsing shape scores for each experiment...")
    # ShapeParser(settings, directories)
    # print("Done.\n")

    
    # print("Imputing missing values for each experiment...")
    # exit_code = subprocess.call("Rscript src/impute_missing_values.r " + config_path, shell = True)
    # if (exit_code > 0):
    #     print("An error occurred while imputing values.")
    #     return()
    # print("Done.\n")

    
    # MachineLearning(directories)


# Entry point
if __name__ == "__main__":
    main()
