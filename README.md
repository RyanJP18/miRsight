- [About](#about)
- [Installation](#installation)
  - [Automatic Setup for Windows/Linux/Mac (Recommended)](#automatic-setup-for-windowslinuxmac-recommended)
    - [Option 1: Docker via CLI](#option-1-docker-via-cli)
    - [Option 2: Docker via GUI](#option-2-docker-via-gui)
  - [Manual Setup (e.g. for HPCs)](#manual-setup-eg-for-hpcs)
- [Usage and Configuration](#usage-and-configuration)
  - [Using Custom Shape and Conservation Data](#using-custom-shape-and-conservation-data)
- [Publications](#publications)
- [Feedback](#feedback)

# About
This is the command line tool used to generate predictions for [miRsight](http://mirsight.info). If you wish to simply use miRsight to produce a set of predictions, it is strongly recommended you use the website directly instead.

If you need to run or configure miRsight locally, miRsight is capable of running on Windows, Linux and Mac on most modern computers via an automated Docker setup. Instructions are also provided for use on HPCs; however, such installations require a degree of manual R and Python package management and are Linux-only due to dependency limitations.

# Installation
 
## Automatic Setup for Windows/Linux/Mac (Recommended) 
1. Download [miRsight](https://github.com/RyanJP18/miRsight/releases)
2. Download and install [Docker](https://www.docker.com/) for your operating system

### Option 1: Docker via CLI
1. Ensure Docker is running-- Windows/Mac: simply run the app, Linux: run `systemctl start docker`
2. Open a terminal and `cd path/to/miRsight`
3. Run `docker build -t mirsight .` to build the image
4. Run `docker run -v output:/app/output -it mirsight` to begin

### Option 2: Docker via GUI
1. Open Docker Desktop and click the `Images` section on the left panel
2. Select `Build`, point to `path/to/miRsight/Dockerfile` and click `Build`
3. Click the `Volumes` section on the left panel
4. Select `Add a File or Folder`, point to `path/to/miRsight/output`
5. In the `Container Path` field, specify `/app/output`
6. Click `Run` to begin

## Manual Setup (e.g. for HPCs)
1. Install the following packages and dependencies:
     - [ViennaRNA 2.5.0](https://github.com/ViennaRNA/ViennaRNA/releases/tag/v2.5.0) (see setup instructions [here](https://github.com/ViennaRNA/ViennaRNA))
     - `R` 4.2.0
       - `jsonlite`
       - `ensembldb`
       - `AnnotationHub`
       - `BSgenome.Hsapiens.NCBI.GRCh38`
       - `dplyr`
       - `stringr`
       - `GenomicScores` (note: only if `use_precompiled_data` is disabled in the config)
     - `python` 3.8
       - `pandas`
       - `pickle`
       - `scikit-learn`
2. Download [miRsight](https://github.com/RyanJP18/miRsight/releases)
3. Open a terminal and `cd path/to/miRsight`
4. run `tar -zxvf precompiled_data.tar.gz` and ensure the `shape` and `output/02-conservation` directories are populated as a result
5. Run `python main.py`

# Usage and Configuration
- Before running miRsight, review the `settings` in `config.json`
- Additionally, edit the `transcript_whitelist.tsv` and `transcript_blacklist.tsv`, and/or `miRNA_whitelist.tsv` and `miRNA_blacklist.tsv`, to avoid unnecessary computation
- Alternatively, provide no filters to generate all miRNA targets in Homo Sapiens; however, this will be slow so you should consider using the website instead
- Output is found in miRsight's `output` folder - if you are only concerned with the final predictions, see `output/10-machine-learning`

## Using Custom Conservation and Shape Data
By default, the `use_precompiled_conservation` and `use_precompiled_shape` flags in `config.json` tell miRsight to use precompiled data. If disabled:

- miRsight will dynamically generate fresh phast100 conservation data against the chosen `ensembl_release` (note: will be slow)
- You can place your own `.shape` output data (from tools like [icSHAPE-pipe](https://github.com/Jun-Lizst/icSHAPE-pipe)) in the `shape` folder to have miRsight use it automatically

# Publications
If you use this tool, please cite: TBD.

# Feedback
If you have any feedback or requests, please submit them via GitHub's issue tracker or contact me at [mail@ryanjp.co.uk](mailto:mail@ryanjp.co.uk).