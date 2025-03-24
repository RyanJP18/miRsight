import csv
import pandas as pd
import numpy as np
import os
from pathlib import Path

class ShapeParser:

    def parse_shape(self, transcript_id, features, shape, shape_source):
        split = features['ensembl_transcript_id_version'].str.split('.').str[0] # get the id without the version
        ts_matches = features.loc[split == transcript_id]

        if (len(ts_matches) != 0): 
            read_len = int(shape[0])
            shape_scores = shape[2:]

            for index, ts in ts_matches.iterrows():
                utr_start = read_len - ts["X3_utr_length"]
                 # match pos is relative to the 6mer, additionally it counts wrong because python goes from 0 whereas R goes from 1. i.e. a match at pos 2 needs to be match_pos - 1 to give it an accessor of 1, or -2 to get the full 8 base seed
                target_start = utr_start + ts["binding_site_pos"] - 2
                target_end = target_start + 8
                sup_start = target_end # 9 (tested, doing a +1 skips a base)
                sup_end = sup_start + 12 # 20
                
                shape_scores_seed_raw = shape_scores[target_start:target_end]
                shape_scores_sup_raw = shape_scores[sup_start:sup_end]
                
                # shape_scores_seed = list(map(float, [score for score in shape_scores_seed_raw if score != "NULL"]))
                # shape_scores_sup = list(map(float, [score for score in shape_scores_sup_raw if score != "NULL"]))

                shape_scores_seed = list(map(float, [raw_score.replace("NULL", "0") for raw_score in shape_scores_seed_raw]))
                shape_scores_sup = list(map(float, [raw_score.replace("NULL", "0") for raw_score in shape_scores_sup_raw]))

                if (len(shape_scores_seed) > 0):
                    features.at[index, shape_source + "_seed"] = np.mean(np.nan_to_num(shape_scores_seed))

                if (len(shape_scores_sup) > 0):
                    features.at[index, shape_source + "_sup"] = np.mean(np.nan_to_num(shape_scores_sup))


    def __init__(self, settings, directories):
        features_files = os.listdir(directories["features_conservation"])

        i = 0
        for features_filename in features_files: 
            i += 1

            if (eval(settings["use_caching"]) and os.path.exists(os.path.join(directories["parsed_shape"], features_filename))):
                print("Shape extraction " + str(i) + "/" + str(len(features_files)) + "... Loaded from cache.")
                continue

            features = pd.read_csv(os.path.join(directories["features_conservation"], features_filename), header = "infer", na_values = "?", sep = "\t")
            features = features[["ensembl_transcript_id_version", "X3_utr_length", "binding_site_pos"]]

            shape_cols = ["ensembl_transcript_id_version"]

            for shape_filename in os.listdir(directories["shape_data"]):
                shape_source = Path(shape_filename).stem.lower()
                features[shape_source + "_seed"] = -1
                features[shape_source + "_sup"] = -1
                # features[shape_source + "_sup_full"] = -1

                shape_cols.append(shape_source + "_seed")
                shape_cols.append(shape_source + "_sup")

                with open(os.path.join(directories["shape_data"], shape_filename)) as shape_file:
                    shapes = csv.reader(shape_file, delimiter="\t")

                    for shape in shapes:
                        self.parse_shape(shape.pop(0), features, shape, shape_source)


            parsed_shape = features[shape_cols]
            parsed_shape.to_csv(os.path.join(directories["parsed_shape"], features_filename), sep = "\t", index = False)

            print("Shape parsing " + str(i) + "/" + str(len(features_files)) + "... Done.")

        