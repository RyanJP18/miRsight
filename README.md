# miRscope
Command line tool used to generate predictions for mirsight.info

# Dependencies (will be simplified in a future release)

- R
- jsonlite
- ensembldb
- AnnotationHub
- BSgenome.Hsapiens.NCBI.GRCh38
- dplyr
- GenomicScores
- stringr

- python
- pandas
- numpy
- pickle
- scikit-learn

# Setup
- Install dependencies
- Unzip "unzip_me.7z" to use pre-compiled SHAPE-seq and conservation data
- Eensure the miRscope/shape and miRscope/output/02-conservation folders are populated
- Run `python main.py'