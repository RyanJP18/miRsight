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
    
    def rank_normalise(self, df, cols):
        """ Convert each cell line's real values to percentiles """
        for c in cols:
            v = pd.to_numeric(df[c], errors="coerce")
            real = v.notna() & (v != -1)
            df[c] = v
            df.loc[real, c] = v[real].rank(pct=True)
        return df

    def combine_weighted(self, df, cols, weights):
        """ Weighted mean of real values """

        w = pd.Series({c: weights[c] for c in cols})
        vals = df[cols].where(df[cols] != -1)
        num = vals.mul(w, axis=1).sum(axis=1, min_count=1)
        den = vals.notna().mul(w, axis=1).sum(axis=1)
        mean = num / den.replace(0, np.nan)
        any_null = (df[cols] == -1).any(axis=1)
        return mean.where(mean.notna(), np.where(any_null, -1.0, np.nan))

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

            feature_headers = next(feature_reader)
            feature_headers.append("shape_seed")
            feature_headers.append("shape_sup")
            feature_headers.append("shape_3")
            feature_headers.append("shape_5")
            feature_headers.append("shape_n_lines")

            features = list(feature_reader)
            shape_df = pd.DataFrame(list(shape_reader))

            groups = [self.shape_seed_cols, self.shape_sup_cols, self.shape_3_cols, self.shape_5_cols]
            
            weights = {c: self.depth[c.split("_")[0]] for g in groups for c in g}
            shape_df = self.rank_normalise(shape_df, [c for g in groups for c in g])

            combined = [self.combine_weighted(shape_df, g, weights) for g in groups]
            real = lambda g: (shape_df[g].notna() & (shape_df[g] != -1)).sum(axis=1)
            n_lines = sum(real(g) for g in groups) / len(groups)

            for i, feature in enumerate(features):
                feature.extend([c.iloc[i] for c in combined])
                feature.append(n_lines.iloc[i])

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
            self.shape_3_cols.append(shape_source + "_3")
            self.shape_5_cols.append(shape_source + "_5")

        with Pool(processes=self.cores) as pool:
            features_files = os.listdir(self.directories["features_conservation"])
            file_count = len(features_files)
            pool.map(self.score_shape, [(features_filename, file_index, file_count) for (file_index, features_filename) in enumerate(features_files)])

    def __init__(self, settings, directories, cores):
        self.settings = settings
        self.directories = directories
        self.cores = int(cores)

        self.use_caching = literal_eval(settings["use_caching"])

        self.depth = {}
        for shape_filename in os.listdir(self.directories["shape_data"]):
            shape_source = Path(shape_filename).stem.lower()
            with open(os.path.join(self.directories["shape_data"], shape_filename)) as f:
                self.depth[shape_source] = sum(1 for _ in f)