"""
Parse shape reactivity values across numerous shape datasets for specific binding site windows.
"""

import csv
import os
from pathlib import Path
from multiprocessing import Pool
from ast import literal_eval
import pandas as pd
import numpy as np


class ShapeParser:

    def parse_shape(self, transcript_id, features, transcript_ids, shape, shape_source):
        """ Parse shape reactivity values for the seed and supplementary portion of a binding from each shape source file """

        ts_matches = features.loc[transcript_ids == transcript_id]
        if len(ts_matches) == 0:
            return

        # the first few cells are metadata, 2+ is shape reactivity scores
        read_len = int(shape[0])
        shape_scores = shape[2:]

        # for each transcript, get shape reactivity scores between specific bases (seed and supplementary portions)
        for index, ts in ts_matches.iterrows():
            utr_start = read_len - ts["X3_utr_length"]
            # note: pos is relative to the 6mer; it also counts wrong because python goes from 0 whereas R goes from 1
            # e.g. a match at pos 2 needs to be match_pos - 1 to give it an accessor of 1, or -2 to get the full 8 base seed
            target_start = utr_start + ts["binding_site_pos"] - 2
            target_end = target_start + 8
            # (note: double tested, doing a +1 skips a base)
            sup_start = target_end
            sup_end = sup_start + 12

            shape_scores_seed_raw = shape_scores[target_start:target_end]
            shape_scores_sup_raw = shape_scores[sup_start:sup_end]

            # shape_scores_seed = list(map(float, [score for score in shape_scores_seed_raw if score != "NULL"]))
            # shape_scores_sup = list(map(float, [score for score in shape_scores_sup_raw if score != "NULL"]))

            shape_scores_seed = list(map(float, [raw_score.replace(
                "NULL", "0") for raw_score in shape_scores_seed_raw]))
            shape_scores_sup = list(map(float, [raw_score.replace(
                "NULL", "0") for raw_score in shape_scores_sup_raw]))

            if len(shape_scores_seed) > 0:
                features.at[index, shape_source +
                            "_seed"] = np.mean(np.nan_to_num(shape_scores_seed))

            if len(shape_scores_sup) > 0:
                features.at[index, shape_source +
                            "_sup"] = np.mean(np.nan_to_num(shape_scores_sup))

    def compute_shape(self, args):
        """ Compute and store shape scores for each shape dataset at each binding site for both the seed and supplementary portions """

        directories, features_filename, file_index, file_count = args

        output_path = os.path.join(
            directories["parsed_shape"], features_filename)

        if self.use_caching and os.path.exists(output_path):
            print(
                f"Shape parse {str(file_index + 1)}/{str(file_count)} - loaded from cache.")
            return

        features = pd.read_csv(os.path.join(
            directories["features_conservation"], features_filename), header="infer", na_values="?", sep="\t")
        features = features[["ensembl_transcript_id_version",
                             "X3_utr_length", "binding_site_pos"]]

        shape_cols = ["ensembl_transcript_id_version"]

        # cycle through each shape file and generate a set of scores for each base
        # note: uses any shape files present in the shape folder and takes and average value between them
        for shape_filename in os.listdir(directories["shape_data"]):
            shape_source = Path(shape_filename).stem.lower()
            features[shape_source + "_seed"] = "NA"
            features[shape_source + "_sup"] = "NA"

            shape_cols.append(shape_source + "_seed")
            shape_cols.append(shape_source + "_sup")

            with open(os.path.join(directories["shape_data"], shape_filename), "r", encoding="utf-8") as shape_file:
                shapes = csv.reader(shape_file, delimiter="\t")

                for shape in shapes:
                    # get the id of the current shape transcript without the version
                    shape_ts_id = shape.pop(0).split(".")[0]
                    transcript_ids = features["ensembl_transcript_id_version"].str.split(
                        ".").str[0]  # get the ids without the version
                    self.parse_shape(shape_ts_id, features,
                                     transcript_ids, shape, shape_source)

        parsed_shape = features[shape_cols]
        parsed_shape.to_csv(output_path, sep="\t", index=False)

        print(f"Shape parse {str(file_index + 1)}/{str(file_count)} - done.")

    def __init__(self, settings, directories, cores):
        self.use_caching = literal_eval(settings["use_caching"])

        with Pool(processes=int(cores)) as pool:
            features_files = os.listdir(directories["features_conservation"])
            file_count = len(features_files)
            pool.map(self.compute_shape, [(directories, features_filename, file_index, file_count) for (
                file_index, features_filename) in enumerate(features_files)])
