"""
Software entry point.
"""

import json
import subprocess
import tarfile
import multiprocessing
import sys
from pathlib import Path
from ast import literal_eval
from src.machine_learning import MachineLearning
from src.conservation_parser import ConservationParser
from src.shape_parser import ShapeParser
from src.shape_scorer import ShapeScorer
from src.rna_folder import RNAFolder


def load_config():
    """ Load and parse the config file to get user settings and directory locations """

    with open(CONFIG_PATH, "r", encoding="utf-8") as config_file:
        config = json.load(config_file)

        dir_paths = list(config["directories"].values())
        for path in dir_paths:
            Path(path).mkdir(parents=True, exist_ok=True)

        max_cores = int(config["settings"]["max_cores"])

        print("The following settings were provided:")
        print(json.dumps(config["settings"], indent=4))
        print(f"See {CONFIG_PATH} to adjust these values.")

        return config["settings"], config["directories"], str(max_cores if max_cores != -1 else multiprocessing.cpu_count() - 1)


def determine_mane_version(ensembl_release):
    """ Load and parse the ensemble:mane lookup file to get ensembl:mane release version mappings """

    with open(ENSEMBL_MANE_LOOKUP_PATH, "r", encoding="utf-8") as mane_ver_file:
        ensembl_mane_lookup = json.load(mane_ver_file)
        return ensembl_mane_lookup[ensembl_release] if ensembl_release in ensembl_mane_lookup else "1.0"


def run_subprocess(command, error_message):
    """ Run a subprocess command and on exception report error details and halt execution """

    try:
        subprocess.run(command, shell=False, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: Command '{e.cmd}' failed with exit code {e.returncode}")
        print(f"Error message: {e.stderr}")
        print(error_message)
        sys.exit(1)


def main(config):
    """ Run each step of the algorithm sequentially to extract and process features in order to ultimately produce predictions """

    settings, directories, cores = config

    print("\n00/11 Downloading annotation data...")
    mane_version = determine_mane_version(settings["ensembl_release"])
    run_subprocess(["sh", "./src/download_annotation_data.sh", settings["use_caching"], directories["preload_data"],
                   settings["ensembl_release"], mane_version], "An error occurred while downloading annotation data.")
    print("00/11 Complete.\n")

    print("01/11 Parsing annotation data and extracting 3' UTR / CDS sequences...")
    run_subprocess(["Rscript", "src/parse_annotation_data.r", CONFIG_PATH], "An error occurred while parsing annotation data.")
    print("01/11 Complete.\n")

    print("02/11 Generating a conservation score cache...")
    if literal_eval(settings["use_precompiled_conservation"]):
        print("Using precompiled data...")
        with tarfile.open(PRECOMPILED_CONSERVATION_PATH, "r:gz") as tar:
            tar.extractall(path=".")
    else:
        print("Using fresh data...")

    run_subprocess(["Rscript", "src/generate_conservation_scores.r", CONFIG_PATH, cores], "An error occurred while generating conservation scores.")
    print("02/11 Complete.\n")

    print("03/11 Locating binding sites for each miRNA...")
    run_subprocess(["Rscript", "src/locate_binding_sites.r", CONFIG_PATH, cores], "An error occurred while locating binding sites.")
    print("03/11 Complete.\n")

    print("04/11 Extracting folding windows for each miRNA...")
    run_subprocess(["Rscript", "src/extract_windows.r", CONFIG_PATH, cores], "An error occurred while extracting folding windows.")
    print("04/11 Complete.\n")

    print("05/11 Folding sequences using ViennaRNA...")
    rna_folder = RNAFolder(settings, directories, cores, 6)

    rna_folder.run_rnafold_batch("windows_rnafold_lr", "folds_rnafold_lr")
    rna_folder.run_rnafold_batch("windows_rnafold_rl", "folds_rnafold_rl")
    rna_folder.run_rnafold_batch("windows_rnafold_ctr", "folds_rnafold_ctr")

    rna_folder.run_rnacofold_batch("windows_rnacofold_full", "folds_rnacofold_full")
    rna_folder.run_rnacofold_batch("windows_rnacofold_seed", "folds_rnacofold_seed")

    rna_folder.run_rnaplfold_batch()
    print("05/11 Complete.\n")

    print("06/11 Extracting features for each miRNA...")
    run_subprocess(["Rscript", "src/extract_features.r", CONFIG_PATH, cores], "An error occurred while extracting features.")
    print("06/11 Complete.\n")

    print("07/11 Parsing conservation scores for each miRNA...")
    conservation_parser = ConservationParser(settings, directories, cores)
    conservation_parser.parse_batch()
    print("07/11 Complete.\n")

    print("08/11 Parsing shape reactivity values for each miRNA...")
    if literal_eval(settings["use_precompiled_shape"]):
        print("Using precompiled data...")
        with tarfile.open(PRECOMPILED_CONSERVATION_PATH, "r:gz") as tar:
            tar.extractall(path=".")
    else:
        print("Using fresh data...")

    shape_parser = ShapeParser(settings, directories, cores)
    shape_parser.parse_batch()
    print("08/11 Complete.\n")

    print("09/11 Producing average shape scores for each miRNA...")
    shape_scorer = ShapeScorer(settings, directories, cores)
    shape_scorer.score_batch()
    print("09/11 Complete.\n")

    print("10/11 Imputing any missing values for each miRNA...")
    run_subprocess(["Rscript", "src/impute_missing_values.r", CONFIG_PATH, cores], "An error occurred while imputing values.")
    print("10/11 Complete.\n")

    print("11/11 Making predictions using machine learning model...")
    machine_learning = MachineLearning(settings, directories, cores)
    machine_learning.bind_model("rf.sav", "scaler.sav")
    machine_learning.predict()
    print("11/11 Complete.\n")


# Entry point
CONFIG_PATH = "config.json"
ENSEMBL_MANE_LOOKUP_PATH = "ensembl_mane_version_lookup.json"
PRECOMPILED_CONSERVATION_PATH = "precompiled_conservation_data.tar.gz"
PRECOMPILED_SHAPE_PATH = "precompiled_shape_data.tar.gz"

if __name__ == "__main__":
    print("\nStarting miRsight...\n")

    main(load_config())

    print("All done!")
    print("See output/11-predictions for the final predictions.\n")
