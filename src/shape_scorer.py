import csv
import pandas as pd
import numpy as np
import os
from pathlib import Path

class ShapeScorer:

    def __init__(self, settings, directories):

        shape_seed_cols = []
        shape_sup_cols = []

        for shape_filename in os.listdir(directories["shape_data"]):
                shape_source = Path(shape_filename).stem.lower()

                shape_seed_cols.append(shape_source + "_seed")
                shape_sup_cols.append(shape_source + "_sup")

        features_files = os.listdir(directories["features_conservation"])
        
        i = 0
        for features_filename in features_files: 
            i += 1

            if (eval(settings["use_caching"]) and os.path.exists(os.path.join(directories["features_cons_shape"], features_filename))):
                print("Shape feature dump " + str(i) + "/" + str(len(features_files)) + "... Loaded from cache.")
                continue

            with open(os.path.join(directories["features_conservation"], features_filename)) as features_file, open(os.path.join(directories["parsed_shape"], features_filename)) as shape_file:
                feature_reader = csv.reader(features_file, delimiter="\t")
                shape_reader = csv.DictReader(shape_file, delimiter="\t")

                feature_headers = next(feature_reader)
                feature_headers.append("shape_seed")
                feature_headers.append("shape_sup")

                features = list(feature_reader)
                parsed_shape = list(shape_reader)

                for feature, shape in zip(features, parsed_shape):
                    total = 0
                    available = 0
                    for col in shape_seed_cols:
                        if float(shape[col]) >= 0:
                            available = available + 1
                            total = total + float(shape[col])
                        
                    if available == 0:
                        feature.append('')
                    else:
                        total = total / available
                        feature.append(total)
                        
                    for col in shape_sup_cols:
                        if float(shape[col]) >= 0:
                            available = available + 1
                            total = total + float(shape[col])
                        
                    if available == 0:
                        feature.append('')
                    else:
                        total = total / available
                        feature.append(total)

                with open(os.path.join(directories["features_cons_shape"], features_filename), "w", newline="") as outfile:
                    writer = csv.writer(outfile, delimiter="\t")
                    writer.writerow(feature_headers)
                    writer.writerows(features)
                    print("Shape feature dump " + str(i) + "/" + str(len(features_files)) + "... Done.")

            