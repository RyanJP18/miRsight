# miRscope
Command line tool used to generate predictions for [miRsight](http://mirsight.info)

# Dependencies 
Note: will be simplified in a future release

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
- Ensure the miRscope/shape and miRscope/output/02-conservation folders are populated as a result
- Run `python main.py`
