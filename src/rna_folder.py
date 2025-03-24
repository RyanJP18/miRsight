import subprocess
import os
from pathlib import Path

class RNAFolder:

    def run_rnafold(self, settings, directories):
        print("RNAfolds 1/3...")
        for experiment_filename in os.listdir(directories["windows_rnafold_lr"]):
            experiment_path = Path(directories["windows_rnafold_lr"], experiment_filename)
            with open(experiment_path) as experiment_file:
                input_file = str(Path(directories["windows_rnafold_lr"], experiment_filename))
                output_file = str(Path(directories["folds_rnafold_lr"], experiment_path.stem))

                if (eval(settings["use_caching"]) and Path(output_file + '.csv').is_file()):
                    print("RNAfold 1/3 " + experiment_path.stem + " already cached, moving to next...")
                else: 
                    with open(f"{output_file}.csv", "w") as out_file:
                        process = subprocess.run(["RNAfold", "--noPS", input_file], stdout=out_file)
                        if (process.returncode != 0):
                            print("An error occurred running RNAfold for " + experiment_filename)
                            continue
        print("Done.\n")

        print("RNAfolds 2/3...")
        for experiment_filename in os.listdir(directories["windows_rnafold_rl"]):
            experiment_path = Path(directories["windows_rnafold_rl"], experiment_filename)
            with open(experiment_path) as experiment_file:
                input_file = str(Path(directories["windows_rnafold_rl"], experiment_filename))
                output_file = str(Path(directories["folds_rnafold_rl"], experiment_path.stem))

                if (eval(settings["use_caching"]) and Path(output_file + '.csv').is_file()):
                    print("RNAfold 2/3 " + experiment_path.stem + " already cached, moving to next...")
                else: 
                    with open(f"{output_file}.csv", "w") as out_file:
                        process = subprocess.run(["RNAfold", "--noPS", input_file], stdout=out_file)
                        if (process.returncode != 0):
                            print("An error occurred running RNAfold for " + experiment_filename)
                            continue
        print("Done.\n")

        print("RNAfolds 3/3...")
        for experiment_filename in os.listdir(directories["windows_rnafold_ctr"]):
            experiment_path = Path(directories["windows_rnafold_ctr"], experiment_filename)
            with open(experiment_path) as experiment_file:
                input_file = str(Path(directories["windows_rnafold_ctr"], experiment_filename))
                output_file = str(Path(directories["folds_rnafold_ctr"], experiment_path.stem))

                if (eval(settings["use_caching"]) and Path(output_file + '.csv').is_file()):
                    print("RNAfold 3/3 " + experiment_path.stem + " already cached, moving to next...")
                else: 
                    with open(f"{output_file}.csv", "w") as out_file:
                        process = subprocess.run(["RNAfold", "--noPS", input_file], stdout=out_file)
                        if (process.returncode != 0):
                            print("An error occurred running RNAfold for " + experiment_filename)
                            continue
        print("Done.\n")

    def run_rnacofold(self, settings, directories):
        print("RNAcofolds 1/2...")
        for experiment_filename in os.listdir(directories["windows_rnacofold_full"]):
            experiment_path = Path(directories["windows_rnacofold_full"], experiment_filename)
            with open(experiment_path) as experiment_file:
                input_file = str(Path(directories["windows_rnacofold_full"], experiment_filename))
                output_file = str(Path(directories["folds_rnacofold_full"], experiment_path.stem))
                
                if (eval(settings["use_caching"]) and Path(output_file + '.csv').is_file()):
                    print("RNAcofold 1/2 " + experiment_path.stem + " already cached, moving to next...")
                else: 
                    with open(f"{output_file}.csv", "w") as out_file:
                        process = subprocess.run(["RNAcofold", "-C", "--noPS", "--output-format=D", input_file], stdout=out_file)
                        if (process.returncode != 0):
                            print("An error occurred running RNAcofold for " + experiment_filename)
                            continue
        print("Done.\n")

        print("RNAcofolds 2/2...")
        for experiment_filename in os.listdir(directories["windows_rnacofold_seed"]):
            experiment_path = Path(directories["windows_rnacofold_seed"], experiment_filename)
            with open(experiment_path) as experiment_file:
                input_file = str(Path(directories["windows_rnacofold_seed"], experiment_filename))
                output_file = str(Path(directories["folds_rnacofold_seed"], experiment_path.stem))

                if (eval(settings["use_caching"]) and Path(output_file + '.csv').is_file()):
                    print("RNAcofold 2/2 " + experiment_path.stem + " already cached, moving to next...")
                else: 
                    with open(f"{output_file}.csv", "w") as out_file:
                        process = subprocess.run(["RNAcofold", "-C", "--noPS", "--output-format=D", input_file], stdout=out_file)
                        if (process.returncode != 0):
                            print("An error occurred running RNAcofold for " + experiment_filename)
                            continue
        print("Done.\n")

    def run_rnaplfold(self, settings, directories):
        print("RNAplfolds...")
        for experiment_filename in os.listdir(directories["windows_rnaplfold"]):
            output_path = Path(directories["folds_rnaplfold"], Path(experiment_filename).stem)
            output_path.mkdir(parents = True, exist_ok = True)
            experiment_path = Path(directories["windows_rnaplfold"], experiment_filename)

            if (eval(settings["use_caching"]) and Path(directories["folds_rnaplfold"], Path(experiment_filename).stem, "sequence_0001_lunp").is_file()):
                print("RNAplfold " + experiment_filename + " already cached, moving to next...")
            else: 
                with open(experiment_path) as experiment_file:
                    input_file = str(Path(directories["windows_rnaplfold"], experiment_filename))
                    exit_code = subprocess.call("RNAplfold -L 40 -W 80 -u 14 --auto-id -o < '" + input_file + "'", shell = True, stderr=subprocess.DEVNULL)
                    # # print(experiment_filename)
                    # with open(input_file, "r") as file:
                    #     content = file.read()

                    #     current_directory = os.getcwd()
                    #     print("Current Directory:", current_directory)

                    #     # process = subprocess.run(["RNAplfold", "-L", "40", "-W", "80", "-u", "14", "--auto-id", "--log-level=3", "-o"], input=content, text=True)
                    #     process = subprocess.run(["RNAplfold", "-L", "40", "-W", "80", "-u", "14", "--auto-id", "-o"], input=content, shell=True, text=True, check=True, capture_output=True)#, stderr=subprocess.DEVNULL)
                        
                        # if (process.returncode != 0):
                    if (exit_code > 0):
                        print("An error occurred running RNAplfold for " + experiment_filename)
                        continue
                    else:
                        if (Path("sequence_0001_lunp").is_file()):
                            subprocess.run(f"mv *_lunp '{output_path}'", shell=True, check=True)
        print("Done.")

    def __init__(self, settings, directories):
        self.run_rnafold(settings, directories)
        self.run_rnacofold(settings, directories)
        self.run_rnaplfold(settings, directories)