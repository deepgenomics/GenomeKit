FROM ubuntu:22.04
COPY genomekit /opt/conda/envs/
ENV PATH=/opt/conda/envs/genomekit/bin:$PATH
