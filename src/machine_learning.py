"""
A utility for preparing a model and necessary features in order to make predictions.
"""

from pathlib import Path
import os
import csv
import pandas as pd

from src.model_trainer import ModelTrainer


class MachineLearning:

    def __init__(self, settings, directories, cores):

        # prepare model
        model = ModelTrainer(cores, directories["features_full_imputed"], directories["features_cons_shape"],
                             directories["machine_learning"], directories["model_data"])

        # load annotations and filters
        annotations = pd.read_csv(
            Path(directories["annotations"], 'annotations.tsv'), sep='\t')
        ensembl_transcript_id_filter = settings['ensembl_transcript_id_filter'].split(
            ',')
        ensembl_gene_id_filter = settings['ensembl_gene_id_filter'].split(',')
        external_gene_id_filter = settings['external_gene_id_filter'].split(
            ',')

        # prepare output files
        output_path_all_pred = Path.joinpath(
            Path(directories["machine_learning"], "all-predictions.tsv"))
        output_path_filtered_pred = Path.joinpath(
            Path(directories["machine_learning"], "filtered-predictions.tsv"))

        if os.path.exists(output_path_all_pred):
            os.remove(output_path_all_pred)

        if os.path.exists(output_path_filtered_pred):
            os.remove(output_path_filtered_pred)

        mirna_files = os.listdir(directories["features_full_imputed"])
        names = [f.split('.')[0] for f in mirna_files]

        i = 0
        # iterate each mirna and make predictions
        for name in names:
            i += 1

            # make predictions
            predictions = model.predict(name)

            # supplement predictions with annotations and reorganise the data to be more useful
            merged_predictions = pd.merge(
                annotations, predictions, on='ensembl_transcript_id_version')
            merged_predictions = merged_predictions[[
                'mirna_id', 'ensembl_transcript_id_version', 'ensembl_gene_id', 'external_gene_id', 'binding_pos', 'seed', 'score']]
            merged_predictions = merged_predictions.sort_values(
                by='score', ascending=False)

            # output unfiltered predictions
            merged_predictions.to_csv(output_path_all_pred, mode='a', sep='\t', index=False, header=(
                i == 1), quoting=csv.QUOTE_NONE)

            # apply any transcript/gene filters supplied in the config
            if settings['ensembl_transcript_id_filter'] != "":
                merged_predictions = merged_predictions[merged_predictions['ensembl_transcript_id_version'].str.split(
                    '.').str[0].isin(ensembl_transcript_id_filter)]
            if settings['ensembl_gene_id_filter'] != "":
                merged_predictions = merged_predictions[merged_predictions['ensembl_gene_id'].isin(
                    ensembl_gene_id_filter)]
            if settings['external_gene_id_filter'] != "":
                merged_predictions = merged_predictions[merged_predictions['external_gene_id'].isin(
                    external_gene_id_filter)]

            # output filtered predictions
            merged_predictions.to_csv(output_path_filtered_pred, mode='a',
                                      sep='\t', index=False, header=(i == 1), quoting=csv.QUOTE_NONE)

            print(f"Predicting targets {i}/{len(names)} - done.")
