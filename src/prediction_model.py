"""
Make predictions for given miRNAs using a trained scikit-learn model.
"""

from pathlib import Path
import pickle
import pandas as pd
import numpy as np
from sklearn.preprocessing import LabelEncoder


class PredictionModel:
    """ A wrapper class that encapsulates a trained scikit-learn model """

    models = None

    def load(self, model_filename):
        """ Load a previous trained scikit learn model from file """

        with open(Path(self.directories["model_data"]).joinpath(model_filename), "rb") as file:
            self.models = pickle.load(file)

    def prep_features(self, data):
        """ Perform any last-minute preprocessing to features ready for prediction, e.g. categorical data encoding """

        dataset = data.copy()

        # drop any unused intermediary features
        dataset = dataset.drop(columns = ["pair_1", "pair_8", "gc_1", "au_8", "au_9", "au_10"], errors="ignore")
        dataset = dataset.drop(columns = ["wobble_pair_1", "wobble_pair_8", "wobble_pair_9", "wobble_pair_10"], errors="ignore")
        dataset = dataset.drop(columns = ["wobble_pair_count_09_20", "wobble_pair_avg_dist_09_20"], errors="ignore")
        dataset = dataset.drop(columns = ["any_pair_count_09_20", "au_count_09_20"], errors="ignore")
        dataset = dataset.drop(columns = ["wobble_pair_count_12_17", "any_pair_count_12_17", "longest_perfect_sequence_12_17",
            "perfect_pair_avg_dist_12_17", "wobble_pair_avg_dist_12_17", "gc_count_12_17", "au_count_12_17", "gc_avg_dist_12_17", "au_avg_dist_12_17",
            "gu_avg_dist_12_17"], errors="ignore") 
        dataset = dataset.drop(columns = ["seed_binding_count",  "longest_perfect_sequence_09_20", "perfect_pair_avg_dist_09_20", 
                                          "seed_gc_ratio", "gc_9", "gc_10", "gc_8", "gc_avg_dist_09_20", "au_avg_dist_09_20", "gu_avg_dist_09_20"], errors="ignore")

        dataset = dataset.drop(columns = ["perfect_pair_count_full", "wobble_pair_count_full", "any_pair_count_full", "longest_perfect_sequence_full", "longest_any_sequence_full", "longest_any_sequence_start_full", "perfect_pair_avg_dist_full", "wobble_pair_avg_dist_full", "any_pair_avg_dist_full",
                                          "gc_count_full", "au_count_full", "gu_count_full", "gc_avg_dist_full", "au_avg_dist_full", "gu_avg_dist_full", "mrna_binding_spread_12_17", "seed_sup_diff_12_17", "mrna_binding_spread_09_20", "seed_sup_diff_09_20", "gc_count_09_20"], errors="ignore")
        dataset = dataset * 1


        seed_map = {"6mer": 0, "7mer-a1": 1, "7mer-m8": 2, "8mer": 3}
        nt_map = {"A": 0, "C": 1, "G": 2, "U": 3}
        dataset['seed_type'] = dataset['seed_binding_type'].map(seed_map)
        dataset['nt_id_at_mirna_1'] = dataset['mirna_1'].map(nt_map)
        dataset['nt_id_at_mirna_8'] = dataset['mirna_8'].map(nt_map)
        dataset = dataset.drop(columns = ["seed_binding_type", "mirna_1", "mirna_8"])#, "mrna_8"])
        
        dataset = dataset.drop(columns = ["binding_site_pos", "au_content_6mer", "site_abundance_7a1", "site_abundance_7m8", "site_abundance_7a1cds", "site_abundance_7m8cds", "au_1",
            "perfect_pair_8", "perfect_pair_count_12_17", "any_pair_avg_dist_12_17", "gu_count_12_17", "longest_any_sequence_09_20", "longest_any_sequence_start_09_20", "seed_sup_diff_full"], errors="ignore")
        dataset = dataset.drop(columns = ["perfect_pair_count_5", "gu_count_5", "longest_any_sequence_5", "longest_any_sequence_start_12_19", "mrna_binding_spread_5", "mrna_binding_spread_12_19",
            "any_pair_avg_dist_12_19", "gu_count_12_19", "perfect_pair_count_12_19"], errors="ignore")

        dataset = dataset.rename(columns={"X3_utr_length": "length_3utr", "cds_length": "length_cds", "site_abundance_7mer": "mts_3utr_7mer", "site_abundance_6mer": "mts_3utr_6mer", "phylo100_seed": "phylop100_seed",
            "site_abundance_6off": "mts_3utr_6off", "shape_n_lines": "shape_contributor_count", "site_abundance_8mer": "mts_3utr_8mer", "site_abundance_6cds": "mts_cds_6mer", "phylo100_3": "phylop100_3",
            "au_content_5_weighted": "au_content_5_wt", "phylo100_5": "phylop100_5", "rnacofold_full_mfe": "rnacofold_mfe_full", "rnacofold_seed_mfe": "rnacofold_mfe_seed", "au_content_3_weighted": "au_content_3_wt", "site_abundance_7cds": "mts_cds_7mer",
            "rel_3utr_pos": "relative_3utr_pos", "phylo100_sup": "phylop100_sup", "mrna_binding_spread_full": "mrna_binding_spread", "best_abundance": "strongest_match", "dist_closest_utr_end": "dist_closest_3utr_end",
            "any_pair_avg_dist_5": "pair_avg_dist_5", "longest_any_sequence_start_5": "longest_seq_start_5", "longest_any_sequence_12_19": "longest_seq_sup", "site_abundance_8cds": "mts_cds_8mer" })

        cols = ["shape_seed", "shape_sup", "shape_3", "shape_5"]
        dataset[cols] = dataset[cols].replace(-1, np.nan)

        print(dataset.columns.tolist())

        return dataset

    def ensemble_predict(self, features):
        """ Rank normalise and produce the mean percentile rank across the ensemble """

        r = np.column_stack([pd.Series(m.predict(features)).rank(pct=True).values for m in self.models])
        return r.mean(axis=1)

    def predict(self, mirna_id):
        """ Produce predictions for the given miRNA """

        raw_test = pd.read_csv(Path.joinpath(Path(self.directories["features_prepared"], mirna_id + ".tsv")), header="infer", na_values="?", sep="\t", index_col=0)
        test = self.prep_features(raw_test)
        test_features = test

        confidence = self.ensemble_predict(test_features)

        # process into a results table
        results = pd.DataFrame(test.index)
        results["mirna_id"] = mirna_id

        raw_test.reset_index(drop=True, inplace=True)
        results.reset_index(drop=True, inplace=True)

        results["score"] = confidence

        # add back in some annotation data
        seeds = raw_test["seed_binding_type"]
        results["seed"] = seeds

        binding_pos = raw_test["binding_site_pos"]
        results["binding_pos"] = binding_pos
        return results.sort_values("score", ascending=False)

    def __init__(self, settings, directories, cores):
        self.settings = settings
        self.directories = directories
        self.cores = int(cores)

        Path(self.directories["machine_learning"]).mkdir(parents=True, exist_ok=True)
