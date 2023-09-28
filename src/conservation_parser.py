import csv
import pandas as pd
import numpy as np
import os
from pathlib import Path

class ConservationParser:

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
            
            start_sup = current_features["binding_site_pos"] + 9
            end_sup = start_sup + 12

            end_3 = start_sup + 30
            end_5 = start_seed - 1
            start_5 = end_5 - 30
            
            conservation_scores_seed_raw = conservation[start_seed:end_seed]
            conservation_scores_sup_raw = conservation[start_sup:end_sup]
            conservation_scores_3_raw = conservation[start_sup:end_3]
            conservation_scores_5_raw = conservation[start_5:end_5]

            # Generally we want to salvage NA values by ignoring them, if more than half of the sequence is NA though log the entire sequence as NA so we can impute the values later
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

    def __init__(self, settings, directories):
        features_files = os.listdir(directories["features"])
        j = 0

        for features_filename in features_files:    
            j += 1

            if (eval(settings["use_caching"]) and os.path.exists(os.path.join(directories["features_conservation"], features_filename))):
                print("Conservation extraction " + str(j) + "/" + str(len(features_files)) + "... Loaded from cache.")
                continue

            features = pd.read_csv(os.path.join(directories["features"], features_filename), header = "infer", na_values = "?", sep = "\t")

            for conservation_filename in os.listdir(directories["conservation"]):
                with open(os.path.join(directories["conservation"], conservation_filename)) as conservation_file:
                    conservations = csv.reader(conservation_file, delimiter=" ")

                    conservation_track_name = Path(conservation_filename).stem
                    features[conservation_track_name + "_seed"] = 0.0
                    features[conservation_track_name + "_sup"] = 0.0
                    
                    i = 0
                    for conservation in conservations:
                        i = self.parse_conservation(i, conservation.pop(0), features, conservation, conservation_track_name)
                        if i == -1:
                            break

                features.to_csv(os.path.join(directories["features_conservation"], features_filename), sep = "\t", index = False)
            
            print("Conservation extraction " + str(j) + "/" + str(len(features_files)) + "... Done.")