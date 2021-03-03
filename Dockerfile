FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /GuestOxenham2021
RUN apt-get update && apt-get install -y --no-install-recommends dirmngr gpg-agent software-properties-common build-essential python3.8 python3-pip python3-setuptools python3-dev
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository -y 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
RUN apt-get -y install r-base
RUN apt-get -y install git
RUN R -e "install.packages(c('merTools', 'dplyr', 'effects', 'ggplot2', 'lme4', 'phia', 'tidyr', 'optimx', 'effects', 'car', 'lmerTest', 'RcppCNPy'), dependencies=TRUE)"
RUN pip3 install numpy scipy pandas scikit-image Cython
RUN pip3 install git+https://github.com/detly/gammatone.git
COPY apcmodels /apcmodels
COPY . /GuestOxenham2021
RUN cd ../apcmodels && python3 setup.py install