FROM tensorflow/tensorflow:latest-gpu
MAINTAINER Cole McKnight <cbmckni@clemson.edu>
MAINTAINER Ben Shealy <btsheal@g.clemson.edu>
MAINTAINER Colin Targonski <ctargon@clemson.edu>

# install package dependencies
RUN apt-get update -qq \
	&& apt-get install -qq -y git python3-pip python-tk

# install python dependencies
RUN pip3 install -q numpy pandas scikit-learn

# install gene-oracle
WORKDIR /opt
RUN git clone https://github.com/SystemsGenetics/gene-oracle.git

ENV GENEORACLE_PATH "/opt/gene-oracle"
ENV PYTHONPATH "$PYTHONPATH:$GENEORACLE_PATH"

# move to home directory
WORKDIR /root
