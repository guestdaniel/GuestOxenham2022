# Start with Ubuntu image
FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive

# Create a directory GuestOxenham2021 and start the container by default in that folder
WORKDIR /GuestOxenham2021

# Install essential tools and Python
RUN apt-get update && apt-get install -y --no-install-recommends dirmngr gpg-agent software-properties-common unzip gcc \
 build-essential python3.8 python3-pip python3-setuptools python3-dev

# Add keys/repos for R 4.0
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository -y 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

# Install R
RUN apt-get -y install r-base

# Install Git
RUN apt-get -y install git

# Install R packages
RUN R -e "install.packages(c('merTools', 'dplyr', 'effects', 'ggplot2', 'lme4', 'phia', 'tidyr', 'optimx', 'effects', 'car', 'lmerTest', 'RcppCNPy'), dependencies=TRUE)"

# Install Python packages
RUN pip3 install numpy scipy pandas scikit-image Cython
RUN pip3 install git+https://github.com/detly/gammatone.git

# Install Verhulst et al. (2018) model
RUN git clone git@github.com:HearingTechnology/Verhulstetal2018Model.git /Verhulstetal2018Model
RUN cd ../Verhulstetal2018Model && gcc -shared -fpic -O3 -ffast-math -o tridiag.so cochlea_utils.c && unzip Poles.zip
RUN export PYTHONPATH="${PYTHONPATH}:/Verhulstetal2018Model"

# Copy apcmodels from disk to image TODO: replace this with a git install of apcmodels
COPY apcmodels /apcmodels

# Copy this folder from disk to image
COPY . /GuestOxenham2021

# Install apcmodels TODO: replace this with a git install of apcmodels
RUN cd ../apcmodels && python3 setup.py install