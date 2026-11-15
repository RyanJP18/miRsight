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

    def parse_shape(self, transcript_id, features, shape, shape_source):
        rows = self._ts_index.get(transcript_id)
        if rows is None:
            return
    
        read_len = int(shape[0])
    
        # parse once per transcript rather than once per window per site.
        scores = np.array([-1.0 if "NULL" in s else float(s) for s in shape[2:]], dtype=float)
        n_scores = scores.size
    
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
    
            for pos in rows:
                utr_start = read_len - features.iat[pos, self._col_utr_len]
                if not (0 <= utr_start < read_len):
                    continue
    
                # match pos is relative to the 6mer, additionally it counts wrong because
                # python goes from 0 whereas R goes from 1. i.e. a match at pos 2 needs to
                # be match_pos - 1 to give it an accessor of 1, or -2 for the full 8 base seed
                target_start = utr_start + features.iat[pos, self._col_bind_pos] - 2
                target_end = target_start + 8
    
                # sup faces miRNA nt 12-19, 5' of the seed match (miRNA 9-11 gap)
                sup_start = target_start - 11
                sup_end = target_start - 3
    
                start_3 = target_end
                end_3 = start_3 + 30
    
                end_5 = target_start
                start_5 = end_5 - 30
    
                if start_5 < 0 or end_3 > n_scores:
                    continue
    
                features.iat[pos, self._col_seed] = np.nanmean(scores[target_start:target_end])
                features.iat[pos, self._col_sup] = np.nanmean(scores[sup_start:sup_end])
                features.iat[pos, self._col_3] = np.nanmean(scores[start_3:end_3])
                features.iat[pos, self._col_5] = np.nanmean(scores[start_5:end_5])

    def parse_source(self, args):
        """ Compute and store shape scores for each shape source by iterating each target row in a features file """

        features_filename, file_index, file_count = args

        output_path = os.path.join(self.directories["parsed_shape"], features_filename)

        if self.use_caching and os.path.exists(output_path):
            print(f"Shape parsing {str(file_index + 1)}/{str(file_count)} - loaded from cache.")
            return

        features = pd.read_csv(os.path.join(self.directories["features_conservation"], features_filename), header="infer", na_values="?", sep="\t")
        features = features[["ensembl_transcript_id_version", "X3_utr_length", "binding_site_pos"]]

        shape_cols = ["ensembl_transcript_id_version"]

        for shape_filename in os.listdir(self.directories["shape_data"]):
            shape_source = Path(shape_filename).stem.lower()
            features[shape_source + "_seed"] = np.nan
            features[shape_source + "_sup"] = np.nan
            features[shape_source + "_3"] = np.nan
            features[shape_source + "_5"] = np.nan
            # features[shape_source + "_sup_full"] = np.nan

            shape_cols.append(shape_source + "_seed")
            shape_cols.append(shape_source + "_sup")
            shape_cols.append(shape_source + "_3")
            shape_cols.append(shape_source + "_5")

            self.prepare_shape_index(features, shape_source)

            with open(os.path.join(self.directories["shape_data"], shape_filename), "r", encoding="utf-8") as shape_file:
                shapes = csv.reader(shape_file, delimiter="\t")

                for shape in shapes:
                    self.parse_shape(shape.pop(0), features, shape, shape_source)

        parsed_shape = features[shape_cols]
        parsed_shape.to_csv(output_path, sep="\t", index=False)

        print(f"Shape parsing {str(file_index + 1)}/{str(file_count)} - done.")

    def parse_batch(self):
        """ Parse shape reactivity values for a batch of features files """

        with Pool(processes=self.cores) as pool:
            features_files = os.listdir(self.directories["features_conservation"])
            file_count = len(features_files)

            pool.map(self.parse_source, [(features_filename, file_index, file_count) for (file_index, features_filename) in enumerate(features_files)])

    def __init__(self, settings, directories, cores):
        self.settings = settings
        self.directories = directories
        self.cores = int(cores)

        self.use_caching = literal_eval(settings["use_caching"])
