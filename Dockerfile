# Start with Ubuntu image
FROM ubuntu:18.04
ENV DEBIAN_FRONTEND=noninteractive

# Create a directory GuestOxenham2021 and start the container by default in that folder
WORKDIR /GuestOxenham2021

# Install essential tools
RUN apt-get update && apt-get install -y --no-install-recommends dirmngr gpg-agent software-properties-common unzip gcc \
    build-essential git libcurl4-openssl-dev

# Add keys/repos for precompiled R packages
RUN apt update -y -qq && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository -y "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    add-apt-repository -y ppa:c2d4u.team/c2d4u4.0+

# Install R
RUN apt-get install -y --no-install-recommends r-base

# Install R packages
#RUN R -e "install.packages(c('merTools', 'dplyr', 'effects', 'ggplot2', 'lme4', 'phia', 'tidyr', 'optimx', 'effects', 'car', 'lmerTest', 'RcppCNPy', 'car'), dependencies=TRUE)"
RUN apt-get install -y -f r-cran-dplyr r-cran-effects r-cran-lme4 r-cran-optimx r-cran-car r-cran-ggplot2
RUN R -e "install.packages(c('merTools', 'phia', 'lmerTest', 'RcppCNPy'), dependencies=TRUE)"

# Install Python
RUN apt-get install -y python3.8 python3-pip python3-setuptools python3-dev

# Install Python packages
RUN pip3 install numpy scipy pandas scikit-image Cython git+https://github.com/detly/gammatone.git
RUN alias python="python3"

# Install Verhulst et al. (2018) model
RUN git clone https://github.com/HearingTechnology/Verhulstetal2018Model /Verhulstetal2018Model
RUN cd ../Verhulstetal2018Model && gcc -shared -fpic -O3 -ffast-math -o tridiag.so cochlea_utils.c
RUN cd ../Verhulstetal2018Model && unzip -v Poles.zip
RUN echo 'export PYTHONPATH="${PYTHONPATH}:/Verhulstetal2018Model"' >> ~/.bashrc

# Copy apcmodels from disk to image TODO: replace this with a git install of apcmodels
COPY apcmodels /apcmodels

# Copy this folder from disk to image
COPY . /GuestOxenham2021

# Install apcmodels TODO: replace this with a git install of apcmodels
RUN cd ../apcmodels && python3 setup.py install