"""
Parse cached conservation scores from generate_conservation_scores at specific target site coordinates by taking the mean of corresponding base values.
"""

import csv
import os
from pathlib import Path
from multiprocessing import Pool
from ast import literal_eval
from collections import namedtuple
import pandas as pd
import numpy as np


class ConservationParser:
    """ A parser which extracts cached generate_conservation_scores output for miRNA target sites """

    def parse_conservation(self, i, transcript_id, features, conservation, conservation_track_name):
        if i >= len(features.index):
            return -1

        current_transcript_id = features.iloc[i]["ensembl_transcript_id_version"]
        if transcript_id != current_transcript_id:

            t_id = float(transcript_id.replace("ENST", ""))
            cur_t_id = float(current_transcript_id.replace("ENST", ""))
            if t_id < cur_t_id:
                return i
            else:
                i += 1
                return self.parse_conservation(i, transcript_id, features, conservation, conservation_track_name)
        
        current_abundance = features.iloc[i]["site_abundance_6mer"]
        for j in range(0, current_abundance):
            current_features = features.iloc[i + j]

            # match pos is relative to the 6mer, additionally it counts wrong because python goes from 0 whereas R goes from 1. i.e. a match at pos 2 needs to be match_pos - 1 to give it an accessor of 1, or -2 to get the full 8 base seed
            start_seed = current_features["binding_site_pos"] - 2
            end_seed = start_seed + 8

            # sup faces miRNA nt 12-19, 5' of the seed match (miRNA 9-11 gap)
            start_sup = start_seed - 11
            end_sup = start_seed - 3

            start_3 = end_seed
            end_3 = start_3 + 30

            end_5 = start_seed
            start_5 = end_5 - 30

            if start_5 < 0 or end_3 > len(conservation):
                continue

            conservation_scores_seed_raw = conservation[start_seed:end_seed]
            conservation_scores_sup_raw = conservation[start_sup:end_sup]
            conservation_scores_3_raw = conservation[start_3:end_3]
            conservation_scores_5_raw = conservation[start_5:end_5]

            # Generally we want to salvage NA values by ignoring them, if more than half of the sequence is NA though log the entire sequence as NA 
            if not "NA" in conservation_scores_seed_raw and len(conservation_scores_seed_raw) > 0:
                conservation_scores_seed = np.array(conservation_scores_seed_raw, dtype=np.float32)
                
                features.at[i + j, conservation_track_name + "_seed"] = np.mean(conservation_scores_seed)
            else:
                features.at[i + j, conservation_track_name + "_seed"] = np.nan


            if not "NA" in conservation_scores_sup_raw and len(conservation_scores_sup_raw) > 0:
                conservation_scores_sup = np.array(conservation_scores_sup_raw, dtype=np.float32)

                features.at[i + j, conservation_track_name + "_sup"] = np.mean(conservation_scores_sup)
            else:
                features.at[i + j, conservation_track_name + "_sup"] = np.nan

            if not "NA" in conservation_scores_3_raw and len(conservation_scores_3_raw) > 0:
                conservation_scores_3 = np.array(conservation_scores_3_raw, dtype=np.float32)
                
                features.at[i + j, conservation_track_name + "_3"] = np.mean(conservation_scores_3)
            else:
                features.at[i + j, conservation_track_name + "_3"] = np.nan


            if not "NA" in conservation_scores_5_raw and len(conservation_scores_5_raw) > 0:
                conservation_scores_5 = np.array(conservation_scores_5_raw, dtype=np.float32)

                features.at[i + j, conservation_track_name + "_5"] = np.mean(conservation_scores_5)
            else:
                features.at[i + j, conservation_track_name + "_5"] = np.nan

        
        return i + j

    def parse_track(self, args):
        """ Compute mean conservation scores for each conservation track by iterating each target row in a features file """

        features_filename, file_index, file_count = args

        output_path = os.path.join(self.directories["features_conservation"], features_filename)

        if self.use_caching and os.path.exists(output_path):
            print(f"Conservation parsing {str(file_index + 1)}/{str(file_count)} - loaded from cache.")
            return

        features = pd.read_csv(os.path.join(self.directories["features"], features_filename), header="infer", na_values="?", sep="\t")

        # cycle through each conservation file and generate a set of scores for each base combination per feature row
        for conservation_filename in os.listdir(self.directories["conservation"]):
            with open(os.path.join(self.directories["conservation"], conservation_filename), "r", encoding="utf-8") as conservation_file:
                conservations = csv.reader(conservation_file, delimiter=" ")

                conservation_track_name = Path(conservation_filename).stem
                features[conservation_track_name + "_seed"] = np.nan
                features[conservation_track_name + "_sup"] = np.nan
                features[conservation_track_name + "_3"] = np.nan
                features[conservation_track_name + "_5"] = np.nan
                
                i = 0
                for conservation in conservations:
                    i = self.parse_conservation(i, conservation.pop(0), features, conservation, conservation_track_name)
                    
                    if i == -1:
                        break

            features.to_csv(output_path, sep="\t", index=False)

        print(f"Conservation parsing {str(file_index + 1)}/{str(file_count)} - done.")

    def parse_batch(self):
        """ Parse conservation scores for a batch of features files """

        with Pool(processes=self.cores) as pool:
            features_files = os.listdir(self.directories["features"])
            file_count = len(features_files)

            pool.map(self.parse_track, [(features_filename, file_index, file_count) for (file_index, features_filename) in enumerate(features_files)])

    def __init__(self, settings, directories, cores):
        self.settings = settings
        self.directories = directories
        self.cores = int(cores)

        self.use_caching = literal_eval(settings["use_caching"])
