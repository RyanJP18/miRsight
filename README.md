- [About](#about)
- [Installation](#installation)
    - [Automatic Setup for Windows/Linux/Mac (Recommended)](#automatic-setup-for-windowslinuxmac-recommended)
        - [Option 1: Docker via CLI](#option-1-docker-via-cli)
        - [Option 2: Docker via GUI](#option-2-docker-via-gui)
    - [Manual Setup (e.g. for HPCs)](#manual-setup-eg-for-hpcs)
- [Configuration](#configuration)
    - [Using Transcript, Gene and miRNA filters](#using-transcript-gene-and-mirna-filters)
    - [Using Custom Conservation and Shape Data](#using-custom-conservation-and-shape-data)
- [Publications](#publications)
- [Feedback](#feedback)

# About
This is the command line tool used to generate mRNA:miRNA target predictions in Homo Sapiens for the [miRsight website](http://mirsight.info). If you wish to simply use miRsight to produce and filter a set of predictions, it is strongly recommended you use the website directly instead.

If you need to run or configure miRsight locally (e.g. to use a specific genome version), miRsight is capable of running on Windows, Linux and Mac on most modern computers via an automated Docker setup. Instructions are also provided for use on HPCs; however, such installations require a degree of manual R and Python package management and are Linux-only due to dependency limitations.

If you use this tool, please cite: TBD.

# Installation
 
### Automatic Setup for Windows/Linux/Mac (Recommended) 
1. Download [miRsight](https://github.com/RyanJP18/miRsight/releases)
2. Download and install [Docker](https://www.docker.com/) for your operating system

##### Option 1: Docker via CLI
1. Ensure Docker is running-- 
   - Windows/Mac: simply run the Docker app
   - Linux: run `systemctl start docker` in a terminal
2. Open a terminal and `cd path/to/miRsight`
3. Run `docker build -t mirsight .` to build the image
4. 
   - Windows: run `docker run -v %cd%/output:/app/output -it mirsight` to begin
   - Linux/Mac: run  `docker run -v $(pwd)/output:/app/output -it mirsight` to begin

##### Option 2: Docker via GUI
1. Open Docker Desktop and click the `Images` section on the left panel
2. Select `Build`, point to `path/to/miRsight/Dockerfile`
3. Click the `Volumes` section on the left panel
4. Select `Add a File or Folder`, point to `path/to/miRsight/output`
5. In the `Container Path` field, specify `/app/output`
6. Click `Run` to begin

### Manual Setup (e.g. for HPCs)
1. Install the following packages and dependencies:
     - [ViennaRNA 2.5.0](https://github.com/ViennaRNA/ViennaRNA/releases/tag/v2.5.0) (see setup instructions [here](https://github.com/ViennaRNA/ViennaRNA))
     - `R 4.2.0`
       - `jsonlite`
       - `ensembldb`
       - `AnnotationHub`
       - `BSgenome.Hsapiens.NCBI.GRCh38`
       - `dplyr`
       - `stringr`
       - `GenomicScores` (note: only if `use_precompiled_conservation` is disabled in the config)
     - `python 3.8`
       - `pandas`
       - `pickle`
       - `scikit-learn`
2. Download [miRsight](https://github.com/RyanJP18/miRsight/releases)
3. Open a terminal and `cd path/to/miRsight`
4. Run `python main.py`

# Configuration
Before running miRsight, review the `settings` in `config.json` and, in particular, consider applying an miRNA filter to speed up computation.

If using the automatic Docker setup, you must rebuild the image after making changes to `config.json`. The image will rebuild much faster than the first time due to caching.

### Using Transcript, Gene and miRNA filters
- One or more comma-separated filters can be provided for:
    - `mirna_id` e.g. `"hsa-miR-129-5p,hsa-miR-30c-5p"`
    - `ensembl_transcript_id` e.g. `"ENST00000000233,ENST00000000412,ENST00000000442"`
    - `ensembl_gene_id` e.g. `"ENSG00000173153,ENSG00000001036"`
    - `external_gene_id` e.g. `"ESRRA,ARF5,SLC7A2,USH1C"`
- A blank filter means a complete whitelist; by default, miRsight will therefore attempt to generate all miRNA targets in Homo Sapiens unless a filter is set
- There is no limit to how stringent or loose these filters can be, but you should consider that looser miRNA filters lead to longer computation times
- Output, both intermediary and final, can found in miRsight's `output` folder - if you are only concerned with the final predictions, see `output/11-target-predictions`

### Using Custom Conservation and Shape Data
By default, the `use_precompiled_conservation` and `use_precompiled_shape` flags in `config.json` tell miRsight to use precompiled data. If disabled:

- miRsight will dynamically generate fresh phylo100 conservation data against the chosen `ensembl_release` (note: will be slow)
- You can place your own `.shape` output data (from tools like [icSHAPE-pipe](https://github.com/Jun-Lizst/icSHAPE-pipe)) in the `shape` folder to have miRsight use it automatically

# Publications
If you use this tool, please cite: TBD.

# Feedback
If you have any feedback or requests, please submit them via GitHub's issue tracker or contact me at [mail@ryanjp.co.uk](mailto:mail@ryanjp.co.uk).
