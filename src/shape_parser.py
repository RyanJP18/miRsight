"""
Parse shape reactivity values across numerous shape datasets for specific binding site windows.
"""

import csv
import os
from pathlib import Path
from multiprocessing import Pool
from ast import literal_eval
from collections import namedtuple
import pandas as pd
import numpy as np


class ShapeParser:
    """ A parser which extracts shape reactivity scores for miRNA target sites """

    ShapeSource = namedtuple("ShapeSource", ["name", "rows"])
    ShapeRow = namedtuple("ShapeRow", ["transcript_id", "read_length", "scores"])

    def _score_target(self, target, shape_row):
        utr_start = shape_row.read_length - target["X3_utr_length"]
        # note: pos is relative to the 6mer; it also counts wrong because python goes from 0 whereas R goes from 1
        # e.g. a match at pos 2 needs to be match_pos - 1 to give it an accessor of 1, or -2 to get the full 8 base seed
        target_start = utr_start + target["binding_site_pos"] - 2
        target_end = target_start + 8
        # (note: double tested, doing a +1 skips a base)
        sup_start = target_end
        sup_end = sup_start + 12

        shape_scores_seed_raw = shape_row.scores[target_start:target_end]
        shape_scores_sup_raw = shape_row.scores[sup_start:sup_end]

        shape_scores_seed = list(map(float, [raw_score.replace("NULL", "0") for raw_score in shape_scores_seed_raw]))
        shape_scores_sup = list(map(float, [raw_score.replace("NULL", "0") for raw_score in shape_scores_sup_raw]))

        seed_score = np.mean(np.nan_to_num(shape_scores_seed)) if len(shape_scores_seed) > 0 else "NA"
        sup_score = np.mean(np.nan_to_num(shape_scores_sup)) if len(shape_scores_sup) > 0 else "NA"

        return seed_score, sup_score

    def _populate_features_with_shape(self, features_with_shape, shape_source):
        """ Parse each shape file for a given features file by walking each row, recursively parsing scores on any targets for which we have data """

        features_with_shape[shape_source.name + "_seed"] = "NA"
        features_with_shape[shape_source.name + "_sup"] = "NA"

        target_transcript_ids = features_with_shape["ensembl_transcript_id_version"].str.split(".").str[0]  # get all target ids without the version number

        for row in shape_source.rows:
            # the first few cells are metadata, 2+ is shape reactivity scores
            shape_transcript_id = row[0].split(".")[0]  # get the id of the current shape transcript without the version number
            read_length = int(row[1])
            shape_scores = row[3:]
            shape_row = self.ShapeRow(shape_transcript_id, read_length, shape_scores)

            target_matches = features_with_shape.loc[target_transcript_ids == shape_row.transcript_id]
            if len(target_matches) == 0:
                continue  # no target matches for this specific shape reactivity row, move to next shape row

            # for each target with a shape reactivity score, get the scores between specific bases (seed and supplementary portions)
            for row_index, target in target_matches.iterrows():
                seed_score, sup_score = self._score_target(target, shape_row)

                features_with_shape.at[row_index, shape_source.name + "_seed"] = seed_score
                features_with_shape.at[row_index, shape_source.name + "_sup"] = sup_score

        return features_with_shape

    def parse_shape(self, args):
        """ Compute and store shape scores for each shape source by iterating each target row in a features file """

        features_filename, file_index, file_count = args

        output_path = os.path.join(self.directories["parsed_shape"], features_filename)

        if self.use_caching and os.path.exists(output_path):
            print(f"Shape parsing {str(file_index + 1)}/{str(file_count)} - loaded from cache.")
            return

        features = pd.read_csv(os.path.join(self.directories["features_conservation"], features_filename), header="infer", na_values="?", sep="\t")
        features_with_shape = features.copy()
        features_with_shape = features_with_shape[["ensembl_transcript_id_version", "X3_utr_length", "binding_site_pos"]]

        shape_cols = ["ensembl_transcript_id_version"]

        # cycle through each shape file and generate a set of scores for each base
        # note: uses any shape files present in the shape folder and takes and average value between them
        for shape_filename in os.listdir(self.directories["shape_data"]):
            with open(os.path.join(self.directories["shape_data"], shape_filename), "r", encoding="utf-8") as shape_file:
                shape_source = self.ShapeSource(Path(shape_filename).stem.lower(), list(csv.reader(shape_file, delimiter="\t")))

                shape_cols.extend([shape_source.name + "_seed", shape_source.name + "_sup"])

                self._populate_features_with_shape(features_with_shape, shape_source)

        parsed_shape = features_with_shape[shape_cols]
        parsed_shape.to_csv(output_path, sep="\t", index=False)

        print(f"Shape parsing {str(file_index + 1)}/{str(file_count)} - done.")

    def parse_batch(self):
        """ Parse shape reactivity values for a batch of features files """

        with Pool(processes=self.cores) as pool:
            features_files = os.listdir(self.directories["features_conservation"])
            file_count = len(features_files)

            pool.map(self.parse_shape, [(features_filename, file_index, file_count) for (file_index, features_filename) in enumerate(features_files)])

    def __init__(self, settings, directories, cores):
        self.settings = settings
        self.directories = directories
        self.cores = int(cores)

        self.use_caching = literal_eval(settings["use_caching"])
