import json
import subprocess
import tarfile
import multiprocessing
import pandas as pd
from pathlib import Path
from src.machine_learning import MachineLearning
from src.conservation_parser import ConservationParser
from src.shape_parser import ShapeParser
from src.shape_scorer import ShapeScorer
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
    with open("ensembl_mane_version_lookup.json") as mane_ver_file:
        ensembl_mane_lookup = json.load(mane_ver_file)
       
        mane_version = "1.0" # default to mane 1.0 if ensembl release specified isn't applicable to try to salvage
        if settings["ensembl_release"] in ensembl_mane_lookup:
            mane_version = ensembl_mane_lookup[settings["ensembl_release"]]
        
        return mane_version

def main():

    print("\n00/11 Downloading annotation data and setting up...")
    config_path = "config.json"
    settings, directories = load_config(config_path)
    mane_version = determine_mane_version(settings)

    max_cores = int(settings['max_cores'])
    cores = str(max_cores if max_cores != -1 else multiprocessing.cpu_count() - 1)
    print("Utilising " + cores + " cores, see config.json to adjust this value.")

    process = subprocess.run(["sh", "./src/download_annotation_data.sh", 
        settings["use_caching"], directories["preload_data"], settings["ensembl_release"], mane_version], shell = False)
    if process.returncode != 0:
        print("An error occurred while downloading annotation data.")
        return()
    print("00/11 Complete.\n")


    print("01/11 Parsing annotation data and extracting 3' UTR / CDS sequences...")
    process = subprocess.run(["Rscript", "src/parse_annotation_data.r", config_path], shell = False)
    if process.returncode != 0:
        print("An error occurred while parsing annotation data.")
        return()
        
    process = subprocess.run(["Rscript", "src/extract_sequences.r", config_path], shell = False)
    if process.returncode != 0:
        print("An error occurred while extracting sequences.")
        return()
    print("01/11 Complete.\n")


    print("02/11 Generating conservation scores...")
    if settings["use_precompiled_conservation"]:
        print("Using precompiled conservation data...")
        with tarfile.open("precompiled_conservation_data.tar.gz", "r:gz") as tar:
            tar.extractall(path=".")

    process = subprocess.run(["Rscript", "src/generate_conservation_scores.r", config_path], shell = False)
    if process.returncode != 0:
        print("An error occurred while generating conservation scores.")
        return()
    print("02/11 Complete.\n")


    print("03/11 Locating binding sites for each miRNA...")
    process = subprocess.run(["Rscript", "src/locate_binding_sites.r", config_path], shell = False)
    if process.returncode != 0:
        print("An error occurred while locating binding sites.")
        return()
    print("03/11 Complete.\n")


    print("04/11 Extracting folding windows for each miRNA...")
    process = subprocess.run(["Rscript", "src/extract_windows.r", config_path], shell = False)
    if process.returncode != 0:
        print("An error occurred while extracting folding windows.")
        return()
    print("04/11 Complete.\n")


    print("05/11 Folding sequences using ViennaRNA...")
    RNAFolder(settings, directories)
    print("05/11 Complete.\n")


    print("06/11 Extracting features for each miRNA...")
    process = subprocess.run(["Rscript", "src/extract_features.r", config_path], shell = False)
    if process.returncode != 0:
        print("An error occurred while extracting features.")
        return()
    print("06/11 Complete.\n")


    print("07/11 Parsing conservation scores for each miRNA...")
    ConservationParser(settings, directories)
    print("07/11 Complete.\n")


    print("08/11 Parsing shape reactivity values for each miRNA...")
    if settings["use_precompiled_shape"]:
        print("Using precompiled shape data...")
        with tarfile.open("precompiled_shape_data.tar.gz", "r:gz") as tar:
            tar.extractall(path=".")

    ShapeParser(settings, directories)
    print("08/11 Complete.\n")


    print("09/11 Producing average shape scores for each miRNA...")
    ShapeScorer(settings, directories)
    print("09/11 Complete.\n")

    
    print("10/11 Imputing any missing values for each miRNA...")
    process = subprocess.run(["Rscript", "src/impute_missing_values.r", config_path], shell = False)
    if process.returncode != 0:
        print("An error occurred while imputing values.")
        return()
    print("10/11 Complete.\n")

    
    print("11/11 Making predictions using machine learning model...")
    MachineLearning(directories)
    print("11/11 Complete.\n")


    print("\nAll done!\n")
    print("See output/11-predictions for the final prediction output.\n")

# Entry point
if __name__ == "__main__":
    main()
