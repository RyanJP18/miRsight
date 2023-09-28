import csv
import pandas as pd
import numpy as np
import os
from pathlib import Path

class ShapeParser:

    def parse_shape(self, transcript_id, features, shape, shape_source):
        ts_matches = features.loc[features['ensembl_transcript_id_version'] == transcript_id]

        if (len(ts_matches) != 0): 
            read_len = int(shape[0])
            shape_scores = shape[2:]

            for index, ts in ts_matches.iterrows():
                utr_start = read_len - ts["X3_utr_length"] 
                #match pos is relative to the 6mer, additionally it counts wrong because python goes from 0 whereas R goes from 1. i.e. a match at pos 2 needs to be match_pos - 1 to give it an accessor of 1, or -2 to get the full 8 base seed
                target_start = utr_start + ts["binding_site_pos"] - 2
                target_end = target_start + 8
                sup_start = target_end # 9
                sup_end = sup_start + 12 # 20
                
                shape_scores_seed_raw = shape_scores[target_start:target_end]
                shape_scores_sup_raw = shape_scores[sup_start:sup_end]

                shape_scores_seed = list(map(float, [raw_score.replace("NULL", "0") for raw_score in shape_scores_seed_raw]))
                shape_scores_sup = list(map(float, [raw_score.replace("NULL", "0") for raw_score in shape_scores_sup_raw]))

                features.at[index, shape_source + "_seed"] = np.mean(np.nan_to_num(shape_scores_seed))
                features.at[index, shape_source + "_sup"] = np.mean(np.nan_to_num(shape_scores_sup))



    def __init__(self, settings, directories):
        features_files = os.listdir(directories["features_conservation"])
        j = 0

        for features_filename in features_files:    
            j += 1

            if (eval(settings["use_caching"]) and os.path.exists(os.path.join(directories["features_cons_shape"], features_filename))):
                print("Shape extraction " + str(j) + "/" + str(len(features_files)) + "... Loaded from cache.")
                continue

            features = pd.read_csv(os.path.join(directories["features_conservation"], features_filename), header = "infer", na_values = "?", sep = "\t")

            for shape_filename in os.listdir(directories["shape_data"]):
                shape_source = Path(shape_filename).stem
                features[shape_source + "_seed"] = -1
                features[shape_source + "_sup"] = -1

                with open(os.path.join(directories["shape_data"], shape_filename)) as shape_file:
                    shapes = csv.reader(shape_file, delimiter="\t")

                    for shape in shapes:
                        self.parse_shape(shape.pop(0), features, shape, shape_source)

            for index, feature in features.iterrows():
                total = 0
                available = 0
                if feature["hela_seed"] >= 0:
                    available = available + 1
                    total = total + feature["hela_seed"]
                
                if feature["hek293_seed"] >= 0:
                    available = available + 1
                    total = total + feature["hek293_seed"]
                
                if feature["h9_seed"] >= 0:
                    available = available + 1
                    total = total + feature["h9_seed"]
                
                if feature["hepg2n_seed"] >= 0:
                    available = available + 1
                    total = total + feature["hepg2n_seed"]
                
                if feature["k562_seed"] >= 0:
                    available = available + 1
                    total = total + feature["k562_seed"]
                
                if available == 0:
                    features.at[index, "shape_seed"] = np.nan
                else:
                    total = total / available
                    features.at[index, "shape_seed"] = total
            
                total = 0
                available = 0
                if feature["hela_sup"] >= 0:
                    available = available + 1
                    total = total + feature["hela_sup"]
                
                if feature["hek293_sup"] >= 0:
                    available = available + 1
                    total = total + feature["hek293_sup"]
                
                if feature["h9_sup"] >= 0:
                    available = available + 1
                    total = total + feature["h9_sup"]
                
                if feature["hepg2n_sup"] >= 0:
                    available = available + 1
                    total = total + feature["hepg2n_sup"]
                
                if feature["k562_sup"] >= 0:
                    available = available + 1
                    total = total + feature["k562_sup"]
                
                if available == 0:
                    features.at[index, "shape_sup"] = np.nan
                else:
                    total = total / available
                    features.at[index, "shape_sup"] = total

            
            features = features.drop(columns = ["hela_seed", "hek293_seed", "h9_seed", "hepg2n_seed", "k562_seed", "hela_sup", "hek293_sup", "h9_sup", "hepg2n_sup", "k562_sup"])

            features.to_csv(os.path.join(directories["features_cons_shape"], features_filename), sep = "\t", index = False)
            print("Shape extraction " + str(j) + "/" + str(len(features_files)) + "... Done.")
            