import csv
import os
import multiprocessing
import pandas as pd
import numpy as np
from pathlib import Path
from multiprocessing import Pool

class ConservationParser:

    def process_score(self, raw_scores):
        # generally we want to salvage NA values by ignoring them, if not possible then log entire sequence as NA so we can impute it later
        if not "NA" in raw_scores and len(raw_scores) > 0:
            return np.mean(np.array(raw_scores, dtype=np.float32))
        else:
            return np.nan

    def parse_transcript(self, i, transcript_id, features, conservation, conservation_track_name):
        if i >= len(features.index):
            return -1 # we have finished looping through each transcript and abundance / mts site

        # if the current id is different to the previous, we are dealing with a new transript
        current_transcript_id = features.iloc[i]["ensembl_transcript_id_version"]
        if transcript_id != current_transcript_id:

            t_id = float(transcript_id.replace("ENST", ""))
            cur_t_id = float(current_transcript_id.replace("ENST", ""))
            if t_id < cur_t_id:
                return i
            else:
                i += 1
                return self.parse_transcript(i, transcript_id, features, conservation, conservation_track_name) # mts site, move index and call recursively

        # mts handling
        current_abundance = features.iloc[i]["site_abundance_6mer"]
        for j in range(0, current_abundance):
            current_features = features.iloc[i + j]

            # match pos is relative to the 6mer, additionally it counts wrong because python goes from 0 whereas R goes from 1
            # i.e. a match at pos 2 needs to be match_pos - 1 to give it an accessor of 1, or -2 to get the full 8 base seed
            start_seed = current_features["binding_site_pos"] - 2 
            end_seed = start_seed + 8
            
            start_sup = current_features["binding_site_pos"] + 9
            end_sup = start_sup + 12

            end_3 = start_sup + 30
            end_5 = start_seed - 1
            start_5 = end_5 - 30
            
            # extract conservation scores between specific base ranges
            features.at[i + j, conservation_track_name + "_seed"] = self.process_score(conservation[start_seed:end_seed])
            features.at[i + j, conservation_track_name + "_sup"] = self.process_score(conservation[start_sup:end_sup])
            features.at[i + j, conservation_track_name + "_3"] = self.process_score(conservation[start_sup:end_3])
            features.at[i + j, conservation_track_name + "_5"] = self.process_score(conservation[start_5:end_5])
            
        return i + j

    def compute_conservation(self, args):
        directories, features_filename, file_index, file_count = args

        output_path = os.path.join(directories["features_conservation"], features_filename)

        if self.use_caching and os.path.exists(output_path):
            print(f"Conservation extraction {str(file_index + 1)}/{str(file_count)} - loaded from cache.")
            return

        features = pd.read_csv(os.path.join(directories["features"], features_filename), header = "infer", na_values = "?", sep = "\t")

        # cycle through each conservation file and generate a set of scores for each base combination per feature row
        # note: the only conservation file used is phylo100 but the mechanism is generic in case more or added in the future (others were used in testing)
        for conservation_filename in os.listdir(directories["conservation"]):
            with open(os.path.join(directories["conservation"], conservation_filename)) as conservation_file:
                conservations = csv.reader(conservation_file, delimiter=" ")

                conservation_track_name = Path(conservation_filename).stem
                features[conservation_track_name + "_seed"] = 0.0
                features[conservation_track_name + "_sup"] = 0.0
                
                row_idx = 0
                for conservation in conservations:
                    row_idx = self.parse_transcript(row_idx, conservation.pop(0), features, conservation, conservation_track_name)
                    if row_idx == -1:
                        break

            features.to_csv(output_path, sep = "\t", index = False)
        
        print(f"Conservation extraction {str(file_index + 1)}/{str(file_count)} - done.")


    def __init__(self, settings, directories):
        self.use_caching = eval(settings["use_caching"])
        max_cores = int(settings['max_cores'])
        cores = max_cores if max_cores != -1 else multiprocessing.cpu_count() - 1

        with Pool(processes=cores) as pool:
            features_files = os.listdir(directories["features"])
            file_count = len(features_files)
            pool.map(self.compute_conservation, [(directories, features_filename, file_index, file_count) for (file_index, features_filename) in enumerate(features_files)])