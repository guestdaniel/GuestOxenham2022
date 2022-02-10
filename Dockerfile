# Start with Ubuntu image
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive

# Create a directory GuestOxenham2021 and start the container by default in that folder
WORKDIR /GuestOxenham2021

# Install essential tools and repositories to get pre-compiled R packages
RUN apt-get update && apt-get install -y --no-install-recommends dirmngr gpg-agent software-properties-common unzip gcc \
    build-essential git libcurl4-openssl-dev openssh-server zlib1g-dev libjpeg-dev && \
    apt update -y -qq && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository -y "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    add-apt-repository -y ppa:c2d4u.team/c2d4u4.0+

# Install R and R packages
RUN apt-get install -y -f --no-install-recommends r-base r-cran-dplyr r-cran-effects r-cran-lme4 r-cran-optimx \
    r-cran-car r-cran-ggplot2 && \
    R -e "install.packages(c('merTools', 'phia', 'lmerTest', 'RcppCNPy'), dependencies=TRUE)"

# Install Python and Python packages
RUN apt-get install -y python3.8 python3-pip python3-setuptools python3-dev && \
    pip3 install numpy==1.19 scipy pandas scikit-image Cython git+https://github.com/detly/gammatone.git && \
    alias python="python3"

# Install Verhulst et al. (2018) model
RUN git clone https://github.com/HearingTechnology/Verhulstetal2018Model /Verhulstetal2018Model && \
    cd ../Verhulstetal2018Model && gcc -shared -fpic -O3 -ffast-math -o tridiag.so cochlea_utils.c && \
    unzip -v Poles.zip && \
    echo 'export PYTHONPATH="${PYTHONPATH}:/Verhulstetal2018Model"' >> ~/.bashrc

# Install apcmodels
RUN git clone https://github.com/guestdaniel/apcmodels /apcmodels && \
    cd ../apcmodels && python3 setup.py install

# Copy this folder from disk to image
COPY . /GuestOxenham2021