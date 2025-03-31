from pathlib import Path
import os
import csv
import multiprocessing

from src.model_trainer import ModelTrainer

class MachineLearning:
    
    def __init__(self, settings, directories):
             
        max_cores = int(settings['max_cores'])
        cores = max_cores if max_cores != -1 else multiprocessing.cpu_count() - 1

        model = ModelTrainer(cores, directories["features_full_imputed"], directories["features_cons_shape"], directories["machine_learning"], directories["model_data"])
        
        experiment_files = os.listdir(directories["features_full_imputed"])
        names = [f.split('.')[0] for f in experiment_files]
        
        i = 0
        for name in names:
            i += 1
            predictions = model.predict(name)
            predictions = predictions.sort_values(by='score')
            predictions.to_csv(Path.joinpath(Path(directories["machine_learning"], "all-predictions.tsv")), mode='a', sep='\t', index=False, quoting=csv.QUOTE_NONE)

            print(f"Predicting targets {i}/{len(names)} - done.")         