"""
Utilise the ViennaRNA suite for RNA folding in order to compute information crucial to binding prediction.
"""

import os
import subprocess
from pathlib import Path
from multiprocessing import Pool
from ast import literal_eval


class RNAFolder:
    """ A utility class for running preset batches of ViennaRNA suite RNA folding tools """

    def run_rnafold_batch(self, input_dir_name, output_dir_name):
        """ Run the RNAfold tool for a given set of input windows to inform accessibility using secondary structure prediction """

        self.current_batch += 1

        input_dir = self.directories[input_dir_name]
        output_dir = self.directories[output_dir_name]

        window_count = len(os.listdir(input_dir))
        i = 0

        for window_filename in os.listdir(input_dir):
            i += 1

            input_path = Path(input_dir, window_filename)
            output_path = Path(output_dir, input_path.stem + ".csv")

            if self.use_caching and output_path.is_file():
                print(
                    f"Folding part {self.current_batch}/{self.batch_count} - {i}/{window_count} - loaded from cache.")
            else:
                exit_code = subprocess.call(
                    f"RNAfold --noPS --jobs={self.cores} {str(input_path)} >> {str(output_path)}", shell=True)
                if exit_code == 0:
                    print(
                        f"Folding part {self.current_batch}/{self.batch_count} - {i}/{window_count} - done.")
                else:
                    print("An error occurred running RNAfold for " + input_path.stem)

    def run_rnacofold_batch(self, input_dir_name, output_dir_name):
        """ Run the RNAcofold tool for a given set of input windows for target binding structure prediction """

        self.current_batch += 1

        input_dir = self.directories[input_dir_name]
        output_dir = self.directories[output_dir_name]

        window_count = len(os.listdir(input_dir))
        i = 0

        for window_filename in os.listdir(input_dir):
            i += 1

            input_path = Path(input_dir, window_filename)
            output_path = Path(output_dir, input_path.stem + ".csv")

            if self.use_caching and output_path.is_file():
                print(
                    f"Folding part {self.current_batch}/{self.batch_count} - {i}/{window_count} - loaded from cache.")
            else:
                exit_code = subprocess.call(
                    f"RNAcofold --jobs={self.cores} -C --noPS --output-format=D {str(input_path)} >> {str(output_path)}", shell=True)
                if exit_code == 0:
                    print(
                        f"Folding part {self.current_batch}/{self.batch_count} - {i}/{window_count} - done.")
                else:
                    print("An error occurred running RNAcofold for " +
                          input_path.stem)

    def run_rnaplfold(self, args):
        """ Run a specific RNAplfold tool instance """

        input_dir, output_dir, window_filename, index, window_count = args

        self.current_batch += 1

        output_path = Path(output_dir, Path(window_filename).stem)

        if self.use_caching and Path(output_path, "sequence_0001_lunp").is_file():
            print(
                f"Folding part {self.current_batch}/{self.batch_count} - {index + 1}/{window_count} - loaded from cache.")
            return

        output_path.mkdir(parents=True, exist_ok=True)
        input_path = Path("..", "..", "..", "..", input_dir, window_filename)

        exit_code = subprocess.call(
            f"RNAplfold -L 40 -W 80 -u 14 --auto-id -o < {str(input_path)}", cwd=output_path, shell=True, stderr=subprocess.DEVNULL)

        if exit_code == 0:
            print(
                f"Folding part {self.current_batch}/{self.batch_count} - {index + 1}/{window_count} - done.")
            subprocess.call("rm -f *_basepairs", cwd=output_path, shell=True)
        else:
            print("An error occurred running RNAplfold for " + input_path.stem)

    def run_rnaplfold_batch(self):
        """ Run the RNAplfold tool for a set of input windows to inform accessibility using secondary structure prediction """

        self.current_batch += 1

        with Pool(processes=self.cores) as pool:
            window_files = os.listdir(self.directories["windows_rnaplfold"])
            window_count = len(window_files)
            
            pool.map(self.run_rnaplfold, [(self.directories["windows_rnaplfold"], self.directories["folds_rnaplfold"],
                     window_filename, index, window_count) for (index, window_filename) in enumerate(window_files)])

    def __init__(self, settings, directories, cores, batch_count):
        self.settings = settings
        self.directories = directories
        self.cores = int(cores)
        self.batch_count = batch_count

        self.current_batch = 0
        self.use_caching = literal_eval(settings["use_caching"])
