# Base Image
FROM python:3.8-slim-buster

# Metadata
LABEL base.image="python:3.8-slim-buster"
LABEL version="1.0"
LABEL software="MemLoci"
LABEL software.version="202001"
LABEL description="an open source software tool to predict subcellular localzation of eukaryotic membrane proteins"
LABEL website="https://github.com/BolognaBiocomp/memloci"
LABEL documentation="https://github.com/BolognaBiocomp/memloci"
LABEL license="GNU GENERAL PUBLIC LICENSE Version 3"
LABEL tags="Proteomics"
LABEL maintainer="Castrense Savojardo <castrense.savojardo2@unibo.it>"

WORKDIR /usr/src/memloci

COPY requirements.txt .

RUN python -m pip install --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt && \
    apt-get -y update && \
    apt-get install -y ncbi-blast+ wget && \
    useradd -m memloci

WORKDIR /seqdb/

RUN wget https://share.biocomp.unibo.it/biocomp/sp2021_01/uniprot_sprot.fasta.gz && \
    gunzip uniprot_sprot.fasta.gz && \
    makeblastdb -in uniprot_sprot.fasta -dbtype prot

WORKDIR /usr/src/memloci

USER memloci

COPY . .

WORKDIR /data/

ENV MEMLOCI_HOME='/usr/src/memloci' PATH=/usr/src/memloci:$PATH

ENTRYPOINT ["/usr/src/memloci/memloci.py"]
