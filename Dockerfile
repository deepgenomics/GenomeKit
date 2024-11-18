FROM ubuntu:22.04
COPY genomekit /opt/conda/envs/genomekit
RUN apt-get update && apt-get install -y wget && rm -rf /var/lib/apt/lists/*
ENV PATH=/opt/conda/envs/genomekit/bin:$PATH
