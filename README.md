# Alignement of phage genome visualisation  
<div style="display" flex; justify-content: space-between; align="center">
  <p>
    <strong>This repository contains scripts used to align and visualise small DNA fragments such as phages genomes.</strong>
  </p>
  <img src="Images/logo_lab.jpg" alt="Lab's logo" width="15%" style="margin-left: 10px;">
</div>

*<div align="center">
    By Ally Champoux, Université de Sherbrooke, 20/12/2024*
</div>

<div align="center">
  
<a href="https://www.gnu.org/software/bash/">![Bash](https://img.shields.io/badge/Shell_script-black?style=for-the-badge&logo=gnubash&logoColor=white)</a>
<a href="https://www.python.org/">![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)</a>

</div>


## Installation 
You have to clone this repository and work in the directory. It contains the following tree: 
```bash
├── LICENSE
├── Test_data
│   └── project
│       ├── data
│       │   ├── blast
│       │   │   └── allPV_vs_allPV.tsv
│       │   └── genomes
│       │       ├── Pv-AC1.gbk
│       │       ├── Pv-AC3.gbk
│       │       ├── Pv-AC4.gbk
│       │       ├── Pv-AM6.gbk
│       │       ├── Pv-AM9.gbk
│       │       ├── Pv-MA11.gbk
│       │       ├── Pv-MA12.gbk
│       │       └── Pv-MA14.gbk
│       ├── genome_comparison_all_PV.svg
│       └── src
│           ├── __init__.py
│           ├── blast_links.py
│           ├── main.py
│           ├── parser.py
│           └── plotter.py
└── README.md

```
You can do so by using the following command: 
```bash
git clone https://github.com/allychamp/Phage_alignement_visualisation.git

```
Note that you should have [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [miniconda](https://docs.anaconda.com/miniconda/) installed to run this pipeline. 
Make sure to activate this environnement to execute the pipeline:
```bash

```
The rest of the dependencies should be installed in the appropriate conda environnement while executing the pipeline.

## Repository content
### `./parser.py`


### `./blast_links.py/`


### `./plotter.py/`

### `./mains.py/`
## Usage
To lunch the script, first complete the config file. Then, go in the directorie containning the script and run the following command:
```bash
python -m src.main
```








Then, the path for the desired ouput directory should be assigned to the `analysis_folder_path` variable just below the `data_folder_path` variable in the snakefile. Once all the paths are set up, make sure to be in the `./phage-genome-analysis-pipeline/` directory and run the pipeline using this command: 
```bash
snakemake --use-conda -j 1  --cores 32 --resources mem_mb=15000
```

The --cores and --resources parameters are set using my computer's resources. Please adapt the command for your computer. In a Linux exploitation system, you can always run :
``` 
free -h
nproc
```
To know excatly how many core and memory are available on your computer. Please use appropriate command for other exploitation systems. Also note that this script is optimised to work with a GPU, it might need adjustments if it is not provided.
## Output
The pipeline will ouput a lot of files. Each sample will have a file looking like this (note that only the important files are represented here. There  a few more that are not shown):
```bash



```


## Tools Used
# Phage_alignement_visualisation
