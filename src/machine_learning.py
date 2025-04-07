"""
Manages machine learning interactions for a trained scikit-learn model and set of data.
"""

from pathlib import Path
import os
import csv
import pandas as pd

from src.prediction_model import PredictionModel


class MachineLearning:
    """ A utility class for managing a trained model and features in order to make predictions """

    model = None
    annotations = None
    mirna_ids = None

    def bind_model(self, model_filename, scaler_filename):
        """ Prepare a trained model and associated scalers for making predictions """

        self.model = PredictionModel(self.settings, self.directories, self.cores)
        self.model.load(model_filename, scaler_filename)

    def predict(self):
        """ Make a set of predictions based on the supplied model and data for each miRNA """

        # wipe any previous predictions
        output_path_all_pred = Path.joinpath(Path(self.directories["machine_learning"], "all-predictions.tsv"))
        output_path_filtered_pred = Path.joinpath(Path(self.directories["machine_learning"], "filtered-predictions.tsv"))

        if os.path.exists(output_path_all_pred):
            os.remove(output_path_all_pred)

        if os.path.exists(output_path_filtered_pred):
            os.remove(output_path_filtered_pred)

        annotations = pd.read_csv(Path(self.directories["annotations"], "annotations.tsv"), sep="\t")

        ensembl_transcript_id_filter = self.settings["ensembl_transcript_id_filter"].split(",")
        ensembl_gene_id_filter = self.settings["ensembl_gene_id_filter"].split(",")
        external_gene_id_filter = self.settings["external_gene_id_filter"].split(",")

        mirna_files = os.listdir(self.directories["features_full_imputed"])
        mirna_ids = [f.split(".")[0] for f in mirna_files]

        # iterate each mirna and make predictions
        for index, mirna_id in enumerate(mirna_ids):
            predictions = self.model.predict(mirna_id)

            # supplement predictions with annotations and reorganise the data to be more useful
            merged_predictions = pd.merge(annotations, predictions, on="ensembl_transcript_id_version")
            merged_predictions = merged_predictions[["mirna_id", "ensembl_transcript_id_version", "ensembl_gene_id", "external_gene_id", "binding_pos", "seed", "score"]]
            merged_predictions = merged_predictions.sort_values(by="score", ascending=False)

            # output unfiltered predictions
            merged_predictions.to_csv(output_path_all_pred, mode="a", sep="\t", index=False, header=(index == 1), quoting=csv.QUOTE_NONE)

            # apply any transcript/gene filters supplied in the config
            if self.settings["ensembl_transcript_id_filter"] != "":
                merged_predictions = merged_predictions[merged_predictions["ensembl_transcript_id_version"].str.split(".").str[0].isin(ensembl_transcript_id_filter)]
            if self.settings["ensembl_gene_id_filter"] != "":
                merged_predictions = merged_predictions[merged_predictions["ensembl_gene_id"].isin(ensembl_gene_id_filter)]
            if self.settings["external_gene_id_filter"] != "":
                merged_predictions = merged_predictions[merged_predictions["external_gene_id"].isin(external_gene_id_filter)]

            # output filtered predictions
            merged_predictions.to_csv(output_path_filtered_pred, mode="a", sep="\t", index=False, header=(index == 1), quoting=csv.QUOTE_NONE)

            print(f"Predicting targets {index + 1}/{len(mirna_ids)} - done.")

    def __init__(self, settings, directories, cores):
        self.settings = settings
        self.directories = directories
        self.cores = cores
