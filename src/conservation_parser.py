"""
Parse cached conservation scores from generate_conservation_scores at specific target site coordinates by taking the mean of corresponding base values.
"""

import csv
import os
from pathlib import Path
from multiprocessing import Pool
from ast import literal_eval
from collections import namedtuple
import pandas as pd
import numpy as np


class ConservationParser:
    """ A utility class for parsing cached generate_conservation_scores output for miRNA target sites """

    ConservationTrack = namedtuple("ConservationTrack", ["name", "rows"])
    Target = namedtuple("Target", ["row_index", "transcript_id"])

    def compute_mean_score(self, raw_scores):
        """ Compute the mean of a collection of scores, or return NA if values are missing """

        return np.nan if "NA" in raw_scores or len(raw_scores) == 0 else np.mean(np.array(raw_scores, dtype=np.float32))

    def score_target(self, index, features, conservation, conservation_track):
        """ Score a specific target by getting the mean across its bases """

        current_features = features.iloc[index]

        # match pos is relative to the 6mer, additionally it counts wrong because python goes from 0 whereas R goes from 1
        # i.e. a match at pos 2 needs to be match_pos - 1 to give it an accessor of 1, or -2 to get the full 8 base seed
        start_seed = current_features["binding_site_pos"] - 2
        end_seed = start_seed + 8

        start_sup = current_features["binding_site_pos"] + 9
        end_sup = start_sup + 12

        end_3 = start_sup + 30
        end_5 = start_seed - 1
        start_5 = end_5 - 30

        # extract conservation scores between specific base ranges
        features.at[index, conservation_track.name + "_seed"] = self.compute_mean_score(conservation[start_seed:end_seed])
        features.at[index, conservation_track.name + "_sup"] = self.compute_mean_score(conservation[start_sup:end_sup])
        features.at[index, conservation_track.name + "_3"] = self.compute_mean_score(conservation[start_sup:end_3])
        features.at[index, conservation_track.name + "_5"] = self.compute_mean_score(conservation[start_5:end_5])

    def parse_row(self, target, features, conservation_row, conservation_track):
        """ Walk a row of the features table, recursively handling any instances of multiple target sites, by parsing and computing mean conservation scores """

        if target.row_index >= len(features.index):  # if our row_index is out of bounds
            return -1  # we have finished looping through each transcript

        # if the current id is different to the previous, we are dealing with a new transript
        current_transcript_id = features.iloc[target.row_index]["ensembl_transcript_id_version"]
        if target.transcript_id != current_transcript_id:

            # compare transcript id numbers to see if we're dealing with multiple target sites or not
            if float(target.transcript_id.replace("ENST", "")) < float(current_transcript_id.replace("ENST", "")):
                return target.row_index  # single target site for this transcript, move on

            # multiple target sites detected, move the index and call recursively
            next_target = self.Target(target.row_index + 1, target.transcript_id)
            return self.parse_row(next_target, features, conservation_row, conservation_track)

        # process each (instance of multiple) target site
        mts_count = features.iloc[target.row_index]["site_abundance_6mer"]
        for i in range(0, mts_count):
            self.score_target(target.row_index + i, features, conservation_row, conservation_track)

        return target.row_index + i

    def parse_conservation_track(self, features, conservation_track):
        """ Parse a conservation track for a given features file """

        features[conservation_track.name + "_seed"] = 0.0
        features[conservation_track.name + "_sup"] = 0.0

        row_idx = 0
        for conservation_row in conservation_track.rows:
            target = self.Target(row_idx, conservation_row.pop(0))

            row_idx = self.parse_row(target, features, conservation_row, conservation_track)
            if row_idx == -1:
                break

    def parse_features_file(self, args):
        """ Compute mean conservation scores for a features file by iterating each of its transcripts with each conservation track """

        features_filename, file_index, file_count = args

        output_path = os.path.join(self.directories["features_conservation"], features_filename)

        if self.use_caching and os.path.exists(output_path):
            print(f"Conservation parsing {str(file_index + 1)}/{str(file_count)} - loaded from cache.")
            return

        features = pd.read_csv(os.path.join(self.directories["features"], features_filename), header="infer", na_values="?", sep="\t")

        # cycle through each conservation file and generate a set of scores for each base combination per feature row
        # note: the only conservation file used is phylo100 but the mechanism is generic (others, such as phast7 and phast100, were also used in testing)
        for conservation_filename in os.listdir(self.directories["conservation"]):
            with open(os.path.join(self.directories["conservation"], conservation_filename), "r", encoding="utf-8") as conservation_file:
                conservation_track = self.ConservationTrack(Path(conservation_filename).stem, list(csv.reader(conservation_file, delimiter=" ")))

                self.parse_conservation_track(features, conservation_track)

                features.to_csv(output_path, sep="\t", index=False)

        print(f"Conservation parsing {str(file_index + 1)}/{str(file_count)} - done.")

    def parse_batch(self):
        """ Parse conservation scores for a batch of files """

        with Pool(processes=self.cores) as pool:
            features_files = os.listdir(self.directories["features"])
            file_count = len(features_files)

            pool.map(self.parse_features_file, [(features_filename, file_index, file_count) for (file_index, features_filename) in enumerate(features_files)])

    def __init__(self, settings, directories, cores):
        self.settings = settings
        self.directories = directories
        self.cores = int(cores)

        self.use_caching = literal_eval(settings["use_caching"])
