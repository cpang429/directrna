FROM nfcore/base
LABEL authors="Tsz Wai Pang" \
      description="Docker image containing all requirements for nf-core/directrna pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
COPY nanoplot-env.yaml /
RUN conda env create -f /nanoplot-env.yaml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-directrna-1.0dev/bin:$PATH
