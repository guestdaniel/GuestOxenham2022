[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4750503.svg)](https://doi.org/10.5281/zenodo.4750503)


# Introduction

This repository contains all the code necessary to completely replicate the figures and analyses found in Guest and Oxenham (2022), in preparation. The codebase is a mixture of Python and R, with Python mostly being used to conduct neural simulations and R mostly being used to analyze behavioral data and create figures. A Docker image is available in which a single script `run.sh` can be run to replicate all of the figures in the paper. This is our recommended solution for running the code in this repository.

The repository includes the code used to collect behavioral data (stimulus generation, stimulus presentation, response collection, etc.) as a single archive. This code is only sparsely documented and relies on MATLAB and the [`afc`](http://medi.uni-oldenburg.de/afc/index.htm) package. The repository does not contain the preprocessed behavioral data but this can be easily downloaded from [Zenodo](https://doi.org/10.5281/zenodo.4750383) and placed in the `data` folder. Code to replicate the analyses and visualizations of the behavioral data and code to replicate the computational models described in the manuscript are available in this respository and are thoroughly documented. The code relies on a number of external libraries, including auditory-specific packages available via our modeling toolbox, [`apcmodels`](https://github.com/guestdaniel/apcmodels).

The code in this repository is licensed under GNU GPLv3 (see `LICENSE.txt`).

# File structure

The file structure of this repository is shown in the file hierarchy below. Essential code required to generate the subfigures in each paper figure is stored in its own folder. A variety of code that was integral to the neural simulations featured in the paper but does not correspond exclusively to one figure is located in the `nofigure` subfolder. Behavioral data analysis scripts are also found in the `nofigure` folder. Finally, a few other scripts and utilities are stored in top-level `.R` or in the `util` folder. Output data files or text files are saved in the same folders as the scripts that generate them, while images/plots are saved in the `plots` folder.

```
.  
├── data                     # Behavioral data goes here!
├── figure1                  # Code for Figure 1
├── ...  
├── figure7                  # Code for Figure 7
│   ├── figure7a.py          # Python code to generate run simulations featured in Figure 7
│   ├── figure7a.R           # R code to plot results from figure7a.py
│   └── ...      
├── figure8                  # Code for Figure 8, similar in structure to /figure7
├── nofigure                 # Code that is needed but is not directly featured in a manuscript figure
│   ├── absolute_thresholds  # Scripts to estimate absolute thresholds at a range of CFs
│   ├── behavioral_data...   # Scripts to perform linear mixed effects modeling on behavioral data
│   ├── tuning_curves        # Scripts to estimate tuning curves and Q10 values at a range of CFs
│   ├── vector_strength      # Scripts to estimate vector strength at a range of CFs
│   └── data_coll...         # .zip file containing code used to collect and preprocess behavioral data
├── plots                    # .png files, exact matches for figures in manuscript
├── supfigure...             # Scripts to generate various supplemental figure components
├── util                     # Miniature package to provide some functions used across repo 
├── config.R                 # Short script to provide constants/configs shared across all R files
├── Dockerfile               # Dockerfile used to generate Docker image
├── LICENSE                  # License file for the code contained in this repository
├── LICENSE_data             # License file for the behavioral data contained in this repository
└── README.md                # This README file
```

# Docker instructions

Docker is our recommended solution for replicating the results from Guest and Oxenham (2022). If you are unfamiliar with Docker, you may want to [orient yourself](https://docs.docker.com/get-started/). The Docker image associated with this repository will allow you to start up a Linux container with Python and R, all required packages/programs, and a copy of this repository. Inside this environment, which is sandboxed from the rest of your system, you can replicate the results of our paper with a single command. When you are done with this process and leave the environment, the environment will clean itself up and stop consuming system resources. The major advantage of using Docker in this way is that you do not have to install Python, R, or any other programs yourself. 

To get started, make sure you have [Docker installed](https://docs.docker.com/get-docker/). Then, follow the instructions below. The instructions below are written for command line interface (such as PowerShell and Terminal) but equivalent commands likely exist in graphical user interface versions of the Docker software.

First, pull the image from our GitHub repository.

```
docker pull docker.pkg.github.com/guestdaniel/guestoxenham2021/guestoxenham2021:1.0.0
```

Next, use the image to create an interactive container.

```
docker run --rm -it guestoxenham2021
```

- `--rm` flag tells Docker to "clean up" after itself and to remove any files generated while running the image after the container is closed
- `-it` tells Docker this is an interactive session 

This container starts with an interactive bash shell located in a copy of the present repository. From there, you can either manually run individual scripts using the `python3` and `Rscript` commands for Python and R, respectively, or you can generate the all the figures in the paper via:
```
bash run.sh
```

However, the figures will be saved out to the container's non-persistent storage and will be destroyed when you exit or end the container. To have permanent copies of the outputs figures saved to your disk, you can link the output `plots` directory inside the container to a preferred output location somewhere on your disk. First, exit the container with the `exit` command, then run the following:

```
docker run --rm -v /home/daniel/GuestOxenham2021/plots:/GuestOxenham2021/plots -it guestoxenham2021
```

- `-v` flag tells Docker to link the `plots` folder on your disk (path to the left of `:`) with the `plots` folder in the container (path to the right of `:`). Obviously, you will need to adjust the path on the  left to point to wherever you have stored your local copy of the repository.

Now, if you call `run.sh` or any of the individual plotting files (e.g., `figure1.R`), whatever is saved in the `plots` folder of the repository will be accessible on your hard drive (outside of the container) in the `output` folder. 

# Data files

## Behavioral data

Behavioral data described in Guest and Oxenham (2022) is stored separately from this code repository on [Zenodo](https://doi.org/10.5281/zenodo.4750383). We recommend that you download the data files (as a `.zip` file) and extract them into the top-level `data` folder in this repository. All the scripts in this repository expect the data to be in this folder. The data is described in a `README.md` file in the Zenodo repository. This step can be performed automatically by the first few lines of `run.sh`. 

# Manual installation

If you do not want to use Docker, you will need Python and R installed on your computer. Instructions for how to install each and configure your environments are provided below. Once you have successfully installed both R and Python and the requisite packages for each, proceed to `Figures` below and read about the code for each figure there. Output figure images will be saved in the `plots` folder as `.png` files.

## Python

A Python 3 interpreter is required to run the simulation code (Figure 6, Figure 6, Figure 6, supplemental figures). We recommend using `pyenv`, `conda`, or another similar tool to install Python 3.6, as well as the packages (with version numbers) listed below:

- `apcmodels` - 0.1.0
- `numpy` - 1.20.21
- `scipy` - 1.6.0
- `pandas` - 1.2.2
- `scikit-image` - 0.18.1
- `Cython` - 0.29.22
- `gammatone`

Presently, `gammatone` can be installed via GitHub as:

```
pip install git+https://github.com/detly/gammatone.git
```

To install `apcmodels`, follow the instructions at https://github.com/guestdaniel/apcmodels.

Once your Python interpreter is configured successfully, set your working directory to your local copy of this repository. Then, run `.py` files as needed. We recommend running the files in the order specified in `run.sh` because some figure files depend on the outputs of earlier scripts.

## R

R is required to generate all the behavioral and modeling figures. The paper figures were generated using R 4.0.3, although in theory any fairly recent version of R should suffice. Below a list of required packages (and the versions used to generate the figures) is provided:

- `merTools` - 0.5.2
- `dplyr` - 1.0.2
- `effects` - 4.2-0
- `ggplot2` - 3.3.2
- `lme4` - 1.1-25
- `phia` - 0.2-1
- `tidyr` - 1.1.2
- `optimx` - 2020-4.2
- `effects` - 4.2-0
- `car` - 3.0-10
- `lmerTest` - 3.1-3
- `RcppCNPy` - 0.2.10

Once your R  is configured successfully, set your working directory to your local copy of this repository. Now, you can run any of the plotting scripts (e.g., `figure1.R`) Some of these scripts depend on the output of the Python scripts. 

# Figures

Specific information about each figure is included below.
Figure generation scripts output figures to the `plots` folder. 
Note that many of the publication figures are generated by combining several subfigures in Inkscape.
Thus, the outputs you should look for in the `plots` folder are generally something like `fig1a.png` for the `.png` file containing Fig 1A. 
The naming scheme is similar but not identical to the naming scheme of the figures in the paper, and some subfigures may be broken up into further subfigures in the `plots` folder (e.g., Figure 5B is composed of `fig5b1.png`, `fig5b2.png`, `fig5b3.png`, and `fig5b4.png`).

### Figure 1

Figure 1 plots example neurograms of harmonic complex tones at low and high frequencies. The entire figure is generated by a single `.py` script.

### Figure 2

Figure 2 plots behavioral results from Experiment 1. The entire figure is generated by a single `.R` script. Corresponding mixed effects models are available in `nofigure/behavioral_data_analysis`. 

### Figure 3

Figure 3 plots behavioral results from Experiment 2. The entire figure is generated by a single `.R` script. Corresponding mixed effects models are available in `nofigure/behavioral_data_analysis`.

### Figure 4

Figure 4 plots excitation patterns for simulated LSR and HSR auditory nerve fibers responding to the ISO and GEOM stimuli from Experiment 1. The simulations and plot are generated by a single `.py` script.

### Figure 5

Figure 5 plots autocorrelograms for the DBL stimuli from Experiment 2. The figure is generated by a single `.py` file.

### Figure 6

Figure 6 plots excitation patterns for the DBL stimuli from Experiment 2. The figure is generated by a single `.py.` file.

### Figure 7

Figure 7 features simulated frequency difference limens (FDLs) derived using ideal observer analysis for three auditory nerve model. Each subfigure has a corresponding `.py` script (to generate the neural simulations) and `.R` scripts (to plot the figure). The `.py` files can take a considerable amount of time and RAM to run, particularly for the sections simulating thresholds for the Verhulst et al. (2018) auditory nerve model.

### Figure 8

Figure 8 features simulated F0 difference limens (F0DLs) derived using ideal observer analysis for three auditory nerve model. Each subfigure has a corresponding `.py` script (to generate the neural simulations) and `.R` scripts (to plot the figure). The `.py` files can take a considerable amount of time and RAM to run, particularly for the sections simulating thresholds for the Verhulst et al. (2018) auditory nerve model.

### Supplemental Figure 1

The supplemental figure features a range of simulation results including vector strength and filter tuning plots for all of the tested auditory nerve models and model responses for the Zilany et al. (2014) auditory nerve model for various types of complex tone stimuli. Each subfigure has a `.py` file to generate it. The first two subfigures rely on simulation results from the `nofigure/vector_strength_curves` and `nofigure/tuning_curves`, respectively. The other subfigures are generated by coresponding `.py` files and can take some time to run. 

### Other supplemental figures

Other supplemental figure generation code is broken up into various scripts with the various `supfigure...` folders. If you have questions about how a particular supplemental figure is generated, please contact the author.