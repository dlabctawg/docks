FROM rocker/verse:3.3.2
MAINTAINER "Brooks Ambrose" brooksambrose@berkeley.edu

## Install igraph
RUN apt-get update && apt-get install -y --no-install-recommends \
	libigraph0 \
	libigraph0-dev

## Install Java and libraries to build rJava
RUN apt-get update && apt-get install -y --no-install-recommends \
	liblzma-dev \
	libbz2-dev \
	libicu-dev \
	default-jdk

## R package prereqs
RUN apt-get update && apt-get install -y --no-install-recommends \
# topicmodels
	libgsl0-dev

### R packages ###

## Tell R where Java is, then test install rJava
RUN R CMD javareconf \
&& Rscript -e 'install.packages(pkgs = c("rJava"), dependencies = TRUE,repos="https://mran.microsoft.com/snapshot/2016-12-25")'

# for jstorr-demo
RUN Rscript -e 'install.packages(pkgs = c("igraph","data.table","formatR","stm","stmBrowser","stmCorrViz","LDAvis","stringdist","qdap","qdapDictionaries","mallet","RTextTools"), dependencies = TRUE,repos="https://mran.microsoft.com/snapshot/2016-12-25")'

RUN Rscript -e 'devtools::install_github("brooksambrose/JSTORr")'

RUN cd ~ && git clone https://github.com/dlabctawg/jstorr-demo

# for chick-or-egg
RUN Rscript -e 'install.packages(pkgs = c("lmtest","forecast","vars"), dependencies = TRUE,repos="https://mran.microsoft.com/snapshot/2016-12-25")'

RUN Rscript -e 'devtools::install_github("seancarmody/ngramr")'

RUN Rscript -e 'install.packages(pkgs = c("magrittr","formatR"), dependencies = TRUE,repos="https://mran.microsoft.com/snapshot/2016-12-25")'

#CMD["\init"]
