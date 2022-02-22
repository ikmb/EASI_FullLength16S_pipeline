FROM nfcore/base
LABEL authors="Eike Matthias Wacker" \
      description="Docker image containing all requirements for EASI FullLength16S pipeline"

COPY environment.yml /
COPY databases /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/ikmb-EASI-pipe-0.1/bin::$PATH
