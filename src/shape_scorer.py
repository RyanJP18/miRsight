import csv
import os
import multiprocessing
import pandas as pd
import numpy as np
from pathlib import Path
from multiprocessing import Pool

class ShapeScorer:

    def compute_average(self, features, parsed_shape):
        for feature, shape in zip(features, parsed_shape):
            total = 0
            available = 0

            # get an average seed reactivity score
            for col in self.shape_seed_cols:
                if shape[col] != 'NA':
                    available = available + 1
                    total = total + float(shape[col])
                
            if available == 0:
                feature.append('NA')
            else:
                total = total / available
                feature.append(total)

            total = 0
            available = 0
                
            # get an average supplementary reactivity score
            for col in self.shape_sup_cols:
                if shape[col] != 'NA':
                    available = available + 1
                    total = total + float(shape[col])
                
            if available == 0:
                feature.append('NA')
            else:
                total = total / available
                feature.append(total)

    def score_shape(self, args):
        directories, features_filename, file_index, file_count = args

        output_path = os.path.join(directories["features_cons_shape"], features_filename)

        if self.use_caching and os.path.exists(output_path):
            print(f"Shape extraction {str(file_index + 1)}/{str(file_count)} - loaded from cache.")
            return

        with open(os.path.join(directories["features_conservation"], features_filename)) as features_file, open(os.path.join(directories["parsed_shape"], features_filename)) as shape_file:
            feature_reader = csv.reader(features_file, delimiter="\t")
            shape_reader = csv.DictReader(shape_file, delimiter="\t")

            # build up a list of new column headers
            feature_headers = next(feature_reader)
            feature_headers.append("shape_seed")
            feature_headers.append("shape_sup")

            # compute the average across each shape column to determine our seed and supplementary average shape reactivity scores
            features = list(feature_reader)
            parsed_shape = list(shape_reader)
            self.compute_average(features, parsed_shape)
        
            with open(output_path, "w", newline="") as outfile:
                writer = csv.writer(outfile, delimiter="\t")
                writer.writerow(feature_headers)
                writer.writerows(features)
                print(f"Shape extraction {str(file_index + 1)}/{str(file_count)} - done.")         

    def __init__(self, settings, directories):
        self.use_caching = eval(settings["use_caching"])
        max_cores = int(settings['max_cores'])
        cores = max_cores if max_cores != -1 else multiprocessing.cpu_count() - 1

        self.shape_seed_cols = []
        self.shape_sup_cols = []

        for shape_filename in os.listdir(directories["shape_data"]):
            shape_source = Path(shape_filename).stem.lower()

            self.shape_seed_cols.append(shape_source + "_seed")
            self.shape_sup_cols.append(shape_source + "_sup")

        with Pool(processes=cores) as pool:
            features_files = os.listdir(directories["features_conservation"])
            file_count = len(features_files)
            pool.map(self.score_shape, [(directories, features_filename, file_index, file_count) for (file_index, features_filename) in enumerate(features_files)])

        