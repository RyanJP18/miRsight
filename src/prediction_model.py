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

    scaler = None
    gs = None
    model = None

    def load(self, model_filename, scaler_filename):
        """ Load a previous trained scikit learn model (and scaler) from file """

        with open(Path(self.directories["model_data"]).joinpath(model_filename), "rb") as file:
            self.gs = pickle.load(file)
            self.model = self.gs.best_estimator_

            if hasattr(self.model, "n_jobs"):
                self.model.set_params(n_jobs=self.cores)

        with open(Path(self.directories["model_data"]).joinpath(scaler_filename), "rb") as file:
            self.scaler = pickle.load(file)

    def prep_features(self, data):
        """ Perform any last-minute preprocessing to features ready for prediction, e.g. categorical data encoding """

        dataset = data.copy()

        # drop unused intermediary features
        dataset = dataset.drop(columns=["seed_binding_count"])
        dataset = dataset.drop(columns=["perfect_pair_count_full", "longest_any_sequence_full", "longest_any_sequence_start_full",
                               "any_pair_avg_dist_full", "gu_count_full", "mrna_binding_spread_12_17", "mrna_binding_spread_09_20"])

        # note: this is now done in imputation step
        # dataset["rel_utr_pos"] = dataset["binding_site_pos"] / dataset["X3_utr_length"]

        dataset = dataset * 1

        # categorical data handling using labelencoder
        labelencoder = LabelEncoder()
        dataset["seed_type"] = labelencoder.fit_transform(dataset["seed_binding_type"])
        dataset["nt_id_at_mirna_1"] = labelencoder.fit_transform(dataset["mirna_1"])
        dataset["nt_id_at_mirna_8"] = labelencoder.fit_transform(dataset["mirna_8"])
        dataset["nt_id_at_mrna_8"] = labelencoder.fit_transform(dataset["mrna_8"])
        dataset = dataset.drop(columns=["seed_binding_type", "mirna_1", "mirna_8", "mrna_8"])

        # note: this is now done in imputation step
        # dataset["X3_utr_length"] = np.log10(dataset["X3_utr_length"])
        # dataset["cds_length"] = np.log10(dataset["cds_length"])
        # dataset["dist_closest_utr_end"] = np.log10(dataset["dist_closest_utr_end"])
        # dataset["binding_site_pos"] = np.log10(dataset["binding_site_pos"])

        dataset = dataset.drop(columns=["binding_site_pos", "site_abundance_7a1", "site_abundance_7m8", "site_abundance_7a1cds", "site_abundance_7m8cds",
                               "perfect_pair_count_12_17", "any_pair_avg_dist_12_17", "gu_count_12_17", "longest_any_sequence_09_20", "longest_any_sequence_start_09_20"])

        return dataset

    def predict(self, mirna_id):
        """ Produce predictions for the given miRNA using the pretrained model and scalers """

        raw_test = pd.read_csv(Path.joinpath(Path(self.directories["features_full_imputed"], mirna_id + ".tsv")), header="infer", na_values="?", sep="\t", index_col=0)
        raw_unimputed_test = pd.read_csv(Path.joinpath(Path(self.directories["features_cons_shape"], mirna_id + ".tsv")), header="infer", na_values="?", sep="\t", index_col=0)

        # manipulate and scale test set features
        test = self.prep_features(raw_test)
        test_features = self.scaler.transform(test)

        # make predictions and get confidence values
        confidence_raw = self.model.predict_proba(test_features)
        conf1d = np.hstack(confidence_raw)
        confidence = conf1d[1::2]

        # process into a results table
        results = pd.DataFrame(test.index)
        results["mirna_id"] = mirna_id

        raw_test.reset_index(drop=True, inplace=True)
        raw_unimputed_test.reset_index(drop=True, inplace=True)
        results.reset_index(drop=True, inplace=True)

        results["score"] = confidence

        # add back in some annotation data
        seeds = raw_test["seed_binding_type"]
        results["seed"] = seeds

        binding_pos = raw_unimputed_test["binding_site_pos"]
        results["binding_pos"] = binding_pos

        # anything above 0.5 confidence we consider to be a prediction
        return results.loc[results["score"] >= 0.5]

    def __init__(self, settings, directories, cores):
        self.settings = settings
        self.directories = directories
        self.cores = int(cores)

        Path(self.directories["machine_learning"]).mkdir(parents=True, exist_ok=True)
