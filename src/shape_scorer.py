"""
Use extracted shape reactivity scores from different sources to produce mean scores for each binding site. 
"""

import csv
import os
from pathlib import Path
from multiprocessing import Pool
from ast import literal_eval


class ShapeScorer:
    """ A parser/scorer which works with intermediary output from shape_parser to compute mean seed and supplementary shape scores across shape sources """

    shape_seed_cols = []
    shape_sup_cols = []

    def _compute_average(self, shape, cols):
        """ Compute the mean shape score by combining each shape source's average """

        # add up the total shape score, noting how many non-NA values were used to build it
        total = 0
        available = 0
        for col in cols:
            if shape[col] != "NA":
                available += 1
                total += float(shape[col])

        # return the average- as NAs are common, try to salvage provided there is at least one value available, otherwise NA
        return "NA" if available == 0 else total / available

    def process_file_pair(self, features, parsed_shape):
        """ Compute the average across each shape column to determine our seed and supplementary average shape reactivity scores for each binding site """

        for feature, shape in zip(features, parsed_shape):
            feature.append(self._compute_average(shape, self.shape_seed_cols))
            feature.append(self._compute_average(shape, self.shape_sup_cols))

        return features

    def score_shape(self, args):
        """ Produce a seed and supplementary shape value for each binding using parsed shape reactivity values from each shape source """

        features_filename, file_index, file_count = args

        output_path = os.path.join(self.directories["features_cons_shape"], features_filename)

        if self.use_caching and os.path.exists(output_path):
            print(f"Shape scoring {str(file_index + 1)}/{str(file_count)} - loaded from cache.")
            return

        with open(os.path.join(self.directories["features_conservation"], features_filename), "r", encoding="utf-8") as features_file, \
                open(os.path.join(self.directories["parsed_shape"], features_filename), "r", encoding="utf-8") as shape_file:

            feature_reader = csv.reader(features_file, delimiter="\t")
            shape_reader = csv.DictReader(shape_file, delimiter="\t")

            # build up a list of new column headers
            feature_headers = next(feature_reader)
            feature_headers.append("shape_seed")
            feature_headers.append("shape_sup")

            features = self.process_file_pair(list(feature_reader), list(shape_reader))

            # store original feature values, plus new shape values in a single table ready for ML
            with open(output_path, "w", encoding="utf-8", newline="") as outfile:
                writer = csv.writer(outfile, delimiter="\t")
                writer.writerow(feature_headers)
                writer.writerows(features)

            print(f"Shape scoring {str(file_index + 1)}/{str(file_count)} - done.")

    def score_batch(self):
        """ Produce mean reactivity scores for a batch of feature files """

        # each shape file located in the main shape directory is considered a new dataset
        for shape_filename in os.listdir(self.directories["shape_data"]):
            shape_source = Path(shape_filename).stem.lower()

            # track a dynamically named (according to the located file) seed and supplementary column for each dataset
            self.shape_seed_cols.append(shape_source + "_seed")
            self.shape_sup_cols.append(shape_source + "_sup")

        with Pool(processes=self.cores) as pool:
            features_files = os.listdir(self.directories["features_conservation"])
            file_count = len(features_files)
            pool.map(self.score_shape, [(features_filename, file_index, file_count) for (
                file_index, features_filename) in enumerate(features_files)])

    def __init__(self, settings, directories, cores):
        self.settings = settings
        self.directories = directories
        self.cores = int(cores)

        self.use_caching = literal_eval(settings["use_caching"])
