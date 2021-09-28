FROM tensorflow/tensorflow:1.15.4-gpu-py3
MAINTAINER Cole McKnight <cbmckni@clemson.edu>
MAINTAINER Ben Shealy <btsheal@clemson.edu>
MAINTAINER Colin Targonski <ctargon@clemson.edu>

# install package dependencies
ENV DEBIAN_FRONTEND=noninteractive 

RUN apt-get update -qq \
    && apt-get install -qq -y git python3-pip

# install python dependencies
RUN pip3 install -q matplotlib numpy pandas scikit-learn seaborn

# install gene-oracle
WORKDIR /opt
RUN git clone -q https://github.com/SystemsGenetics/gene-oracle.git

ENV GENEORACLE_PATH "/opt/gene-oracle"
ENV PYTHONPATH "${PYTHONPATH}:${GENEORACLE_PATH}"

# initialize default work directory
WORKDIR /workspace
