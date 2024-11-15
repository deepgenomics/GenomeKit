FROM ubuntu:22.04
COPY genomekit /opt/conda/envs/genomekit
ENV PATH=/opt/conda/envs/genomekit/bin:$PATH
