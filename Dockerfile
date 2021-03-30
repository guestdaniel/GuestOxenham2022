# Start with Ubuntu image
FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive

# Create a directory GuestOxenham2021 and start the container by default in that folder
WORKDIR /GuestOxenham2021

# Install essential tools
RUN apt-get update && apt-get install -y --no-install-recommends dirmngr gpg-agent software-properties-common unzip gcc \
 build-essential git

# Add keys/repos for R 4.0
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository -y 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

# Install R
RUN apt-get -y install r-base

# Install R packages
RUN R -e "install.packages(c('merTools', 'dplyr', 'effects', 'ggplot2', 'lme4', 'phia', 'tidyr', 'optimx', 'effects', 'car', 'lmerTest', 'RcppCNPy'), dependencies=TRUE)"

# Install Python
RUN apt-get install -y python3.8 python3-pip python3-setuptools python3-dev

# Install Python packages
RUN pip3 install numpy scipy pandas scikit-image Cython
RUN pip3 install git+https://github.com/detly/gammatone.git
RUN alias python="python3"

# Install Verhulst et al. (2018) model
RUN git clone https://github.com/HearingTechnology/Verhulstetal2018Model /Verhulstetal2018Model
RUN cd ../Verhulstetal2018Model && gcc -shared -fpic -O3 -ffast-math -o tridiag.so cochlea_utils.c
RUN cd ../Verhulstetal2018Model && unzip -v Poles.zip
RUN export PYTHONPATH="${PYTHONPATH}:/Verhulstetal2018Model"

# Copy apcmodels from disk to image TODO: replace this with a git install of apcmodels
COPY apcmodels /apcmodels

# Copy this folder from disk to image
COPY . /GuestOxenham2021

# Install apcmodels TODO: replace this with a git install of apcmodels
RUN cd ../apcmodels && python3 setup.py install