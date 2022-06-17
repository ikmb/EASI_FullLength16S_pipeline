FROM nfcore/base
LABEL authors="Eike Matthias Wacker" \
      description="Docker image containing all requirements for EASI FullLength16S pipeline"

COPY environment.yml /
#COPY databases /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/ikmb-EASI-pipe-0.1/bin::$PATH

RUN wget https://zenodo.org/record/4310151/files/rdp_species_assignment_18.fa.gz
RUN wget https://zenodo.org/record/4310151/files/rdp_train_set_18.fa.gz
RUN wget https://zenodo.org/record/4735821/files/GTDB_bac120_arc122_ssu_r202_Genus.fa.gz
RUN wget https://zenodo.org/record/4735821/files/GTDB_bac120_arc122_ssu_r202_Species.fa.gz