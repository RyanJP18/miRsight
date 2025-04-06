"""
Produce representative shape scores for the seed and supplementary portion of each binding site. 
"""

import csv
import os
from pathlib import Path
from multiprocessing import Pool
from ast import literal_eval


class ShapeScorer:

    def compute_average(self, shape, cols):
        """ Compute the mean shape score by combining each shape source's average """

        total = 0
        available = 0

        # get the total reactivity score for the provided columns
        for col in cols:
            if shape[col] != 'NA':
                available = available + 1
                total = total + float(shape[col])

        # return the average- as NAs are common, try to salvage provided there is at least one value available, otherwise NA
        return 'NA' if available == 0 else total / available

    def score_shape(self, args):
        """ Produce a seed and supplementary shape value for each binding using parsed shape reactivity values from each shape source """

        directories, features_filename, file_index, file_count = args

        output_path = os.path.join(
            directories["features_cons_shape"], features_filename)

        if self.use_caching and os.path.exists(output_path):
            print(
                f"Shape extraction {str(file_index + 1)}/{str(file_count)} - loaded from cache.")
            return

        with open(os.path.join(directories["features_conservation"], features_filename), "r", encoding="utf-8") as features_file, \
                open(os.path.join(directories["parsed_shape"], features_filename), "r", encoding="utf-8") as shape_file:

            feature_reader = csv.reader(features_file, delimiter="\t")
            shape_reader = csv.DictReader(shape_file, delimiter="\t")

            # build up a list of new column headers
            feature_headers = next(feature_reader)
            feature_headers.append("shape_seed")
            feature_headers.append("shape_sup")

            # compute the average across each shape column to determine our seed and supplementary average shape reactivity scores
            features = list(feature_reader)
            parsed_shape = list(shape_reader)
            for feature, shape in zip(features, parsed_shape):
                feature.append(self.compute_average(
                    shape, self.shape_seed_cols))
                feature.append(self.compute_average(
                    shape, self.shape_sup_cols))

            # store original feature values, plus new shape values in a single table ready for ML
            with open(output_path, "w", encoding="utf-8", newline="") as outfile:
                writer = csv.writer(outfile, delimiter="\t")
                writer.writerow(feature_headers)
                writer.writerows(features)
                print(
                    f"Shape extraction {str(file_index + 1)}/{str(file_count)} - done.")

    def __init__(self, settings, directories, cores):
        self.use_caching = literal_eval(settings["use_caching"])

        self.shape_seed_cols = []
        self.shape_sup_cols = []

        for shape_filename in os.listdir(directories["shape_data"]):
            shape_source = Path(shape_filename).stem.lower()

            self.shape_seed_cols.append(shape_source + "_seed")
            self.shape_sup_cols.append(shape_source + "_sup")

        with Pool(processes=int(cores)) as pool:
            features_files = os.listdir(directories["features_conservation"])
            file_count = len(features_files)
            pool.map(self.score_shape, [(directories, features_filename, file_index, file_count) for (
                file_index, features_filename) in enumerate(features_files)])
