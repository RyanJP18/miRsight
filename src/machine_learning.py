from pathlib import Path
import os
import csv

from src.model_trainer import ModelTrainer

class MachineLearning:
    
    def __init__(self, directories):
             
        model = ModelTrainer(directories["features_full_imputed"], directories["features_cons_shape"], directories["machine_learning"], directories["model_data"])
        
        experiment_files = os.listdir(directories["features_full_imputed"])
        names = [f.split('.')[0] for f in experiment_files]
        for name in names:
            predictions = model.predict(name)
            predictions = predictions.sort_values(by='score')
            predictions.to_csv(Path.joinpath(Path(directories["machine_learning"], "all-predictions.tsv")), mode='a', sep='\t', index=False, quoting=csv.QUOTE_NONE)