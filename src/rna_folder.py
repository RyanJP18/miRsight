import os
import subprocess
import multiprocessing
import logging
from pathlib import Path
from multiprocessing import Pool

class RNAFolder:

    def run_rnafold(self, input_dir, output_dir, step, use_caching, cores):
        window_count = len(os.listdir(input_dir))
        i = 0

        for window_filename in os.listdir(input_dir):
            i += 1

            input_path = Path(input_dir, window_filename)
            output_path = Path(output_dir, input_path.stem + '.csv')

            if use_caching and output_path.is_file():
                print(f"Folding part {step}/6 - {i}/{window_count} - loaded from cache.")
            else: 
                exit_code = subprocess.call(f"RNAfold --noPS --jobs={cores} {str(input_path)} >> {str(output_path)}", shell = True)
                if exit_code == 0:
                    print(f"Folding part {step}/6 - {i}/{window_count} - done.")
                else:
                    print("An error occurred running RNAfold for " + input_path.stem)

    def run_rnacofold(self, input_dir, output_dir, step, use_caching, cores):
        window_count = len(os.listdir(input_dir))
        i = 0

        for window_filename in os.listdir(input_dir):
            i += 1

            input_path = Path(input_dir, window_filename)
            output_path = Path(output_dir, input_path.stem + '.csv')
                
            if use_caching and output_path.is_file():
                print(f"Folding part {step}/6 - {i}/{window_count} - loaded from cache.")
            else: 
                exit_code = subprocess.call(f"RNAcofold --jobs={cores} -C --noPS --output-format=D {str(input_path)} >> {str(output_path)}", shell = True)
                if exit_code == 0:
                    print(f"Folding part {step}/6 - {i}/{window_count} - done.")
                else:
                    print("An error occurred running RNAcofold for " + input_path.stem)
        
    def run_rnaplfold(self, args):
        input_dir, output_dir, window_filename, index, window_count, use_caching = args

        output_path = Path(output_dir, Path(window_filename).stem)

        if use_caching and Path(output_path, "sequence_0001_lunp").is_file():
            print(f"Folding part 6/6 - {index + 1}/{window_count} - loaded from cache.")
            return
            
        output_path.mkdir(parents = True, exist_ok = True)
        input_path = Path("..", "..", "..", "..", input_dir, window_filename)

        exit_code = subprocess.call(f"RNAplfold -L 40 -W 80 -u 14 --auto-id -o < {str(input_path)}", cwd=output_path, shell = True, stderr=subprocess.DEVNULL)
    
        if exit_code == 0:
            print(f"Folding part 6/6 - {index + 1}/{window_count} - done.")
            subprocess.call("rm -f *_basepairs", cwd=output_path, shell=True)
        else:
            print("An error occurred running RNAplfold for " + input_path.stem)

    def __init__(self, settings, directories):
        use_caching = eval(settings["use_caching"])
        max_cores = int(settings['max_cores'])
        cores = max_cores if max_cores != -1 else multiprocessing.cpu_count() - 1

        self.run_rnafold(directories["windows_rnafold_lr"], directories["folds_rnafold_lr"], 1, use_caching, cores)
        self.run_rnafold(directories["windows_rnafold_rl"], directories["folds_rnafold_rl"], 2, use_caching, cores)
        self.run_rnafold(directories["windows_rnafold_ctr"], directories["folds_rnafold_ctr"], 3, use_caching, cores)
        self.run_rnacofold(directories["windows_rnacofold_full"], directories["folds_rnacofold_full"], 4, use_caching, cores)
        self.run_rnacofold(directories["windows_rnacofold_seed"], directories["folds_rnacofold_seed"], 5, use_caching, cores)
        
        with Pool(processes=cores) as pool:
            window_files = os.listdir(directories["windows_rnaplfold"])
            window_count = len(window_files)
            pool.map(self.run_rnaplfold, [(directories["windows_rnaplfold"], directories["folds_rnaplfold"], window_filename, index, window_count, use_caching) for (index, window_filename) in enumerate(window_files)])