# Base Image
FROM continuumio/miniconda3:4.8

# Metadata
LABEL base.image="continuumio/miniconda3:4.8"
LABEL version="3.0"
LABEL software="TPpred3"
LABEL software.version="2015012"
LABEL description="an open source software tool to predict organelle-targeting peptides in proteins"
LABEL website="https://tppred3.biocomp.unibo.it"
LABEL documentation="https://tppred3.biocomp.unibo.it"
LABEL license="GNU GENERAL PUBLIC LICENSE Version 3"
LABEL tags="Proteomics"
LABEL maintainer="Castrense Savojardo <castrense.savojardo2@unibo.it>"

ENV PYTHONDONTWRITEBYTECODE=true TPPRED_ROOT=/usr/src/tppred3 PATH=/usr/src/tppred3:$PATH

WORKDIR /usr/src/tppred3

COPY . .

WORKDIR /data/

RUN conda update -n base conda && \
   conda install --yes nomkl meme emboss -c bioconda && \
   conda install --yes nomkl libsvm -c conda-forge && \
   conda install --yes nomkl libiconv biopython && \
   conda clean -afy \
   && find /opt/conda/ -follow -type f -name '*.a' -delete \
   && find /opt/conda/ -follow -type f -name '*.pyc' -delete \
   && find /opt/conda/ -follow -type f -name '*.js.map' -delete

ENTRYPOINT ["/usr/src/tppred3/tppred3.py"]
