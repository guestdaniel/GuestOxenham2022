FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /GuestOxenham2021
RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base python3.8 python3-pip python3-setuptools python3-dev git
RUN R -e "install.packages(c('merTools', 'dplyr', 'effects', 'ggplot2', 'lme4', 'phia', 'tidyr', 'optimx', 'effects', 'car', 'lmerTest', 'RcppCNPy'), dependencies=TRUE)"
RUN pip3 install numpy scipy pandas scikit-image Cython
RUN pip3 install git+https://github.com/detly/gammatone.git
COPY apcmodels /apcmodels
COPY . /GuestOxenham2021
RUN cd ../apcmodels && python3 setup.py install