import os
import subprocess
import multiprocessing
import logging
from pathlib import Path
from multiprocessing import Pool

class RNAFolder:

    def run_rnafold(self, input_dir, output_dir, use_caching, cores):
        for window_filename in os.listdir(input_dir):
            input_path = Path(input_dir, window_filename)
            output_path = Path(output_dir, input_path.stem + '.csv')

            if use_caching and output_path.is_file():
                print(input_path.stem + " already cached, moving to next...")
            else: 
                exit_code = subprocess.call(f"RNAfold --noPS --jobs={cores} {str(input_path)} >> {str(output_path)}", shell = True)
                if exit_code == 0:
                    print(input_path.stem + " folding complete.")
                else:
                    print("An error occurred running RNAfold for " + input_path.stem)

    def run_rnacofold(self, input_dir, output_dir, use_caching, cores):
        for window_filename in os.listdir(input_dir):
            input_path = Path(input_dir, window_filename)
            output_path = Path(output_dir, input_path.stem + '.csv')
                
            if use_caching and output_path.is_file():
                print(input_path.stem + " already cached, moving to next...")
            else: 
                exit_code = subprocess.call(f"RNAcofold --jobs={cores} -C --noPS --output-format=D {str(input_path)} >> {str(output_path)}", shell = True)
                if exit_code == 0:
                    print(input_path.stem + " folding complete.")
                else:
                    print("An error occurred running RNAcofold for " + input_path.stem)
        
    def run_rnaplfold(self, args):
        input_dir, output_dir, window_filename, use_caching = args

        output_path = Path(output_dir, Path(window_filename).stem)
        output_path.mkdir(parents = True, exist_ok = True)
        os.chdir(output_path)

        input_path = Path("..", "..", "..", "..", input_dir, window_filename)

        if use_caching and Path(output_path, "sequence_0001_lunp").is_file():
            print(input_path.stem + " already cached, moving to next...", flush=True)
        else: 
            exit_code = subprocess.call(f"RNAplfold -L 40 -W 80 -u 14 --auto-id -o < {str(input_path)}", shell = True, stderr=subprocess.DEVNULL)
            if exit_code == 0:
                print(input_path.stem + " folding complete.")
                if Path("sequence_0001_basepairs").is_file():
                    subprocess.run(f"rm -f *_basepairs", shell=True, check=True)
            else:
                print("An error occurred running RNAplfold for " + input_path.stem, flush=True)
        

    def __init__(self, settings, directories):
        use_caching = eval(settings["use_caching"])
        max_cores = int(settings['max_cores'])
        cores = max_cores if max_cores != -1 else multiprocessing.cpu_count() - 1

        print("RNAfolds 1/6...")
        self.run_rnafold(directories["windows_rnafold_lr"], directories["folds_rnafold_lr"], use_caching, cores)
        print("Done.")

        print("RNAfolds 2/6...")
        self.run_rnafold(directories["windows_rnafold_rl"], directories["folds_rnafold_rl"], use_caching, cores)
        print("Done.")

        print("RNAfolds 3/6...")
        self.run_rnafold(directories["windows_rnafold_ctr"], directories["folds_rnafold_ctr"], use_caching, cores)
        print("Done.")


        print("RNAcofolds 4/6...")
        self.run_rnacofold(directories["windows_rnacofold_full"], directories["folds_rnacofold_full"], use_caching, cores)
        print("Done.")

        print("RNAcofolds 5/6...")
        self.run_rnacofold(directories["windows_rnacofold_seed"], directories["folds_rnacofold_seed"], use_caching, cores)
        print("Done.")
        
        print("RNAplfolds 6/6...")
        with Pool(processes=cores) as pool:
            pool.map(self.run_rnaplfold, [(directories["windows_rnaplfold"], directories["folds_rnaplfold"], window_filename, use_caching) for window_filename in os.listdir(directories["windows_rnaplfold"])])
        print("Done.")