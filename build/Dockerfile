FROM rocker/verse:3.3.2
MAINTAINER "Brooks Ambrose" brooksambrose@berkeley.edu

## Add repository for java 8
RUN apt-get install software-properties-common \
	&& echo "deb http://ppa.launchpad.net/webupd8team/java/ubuntu xenial main" | tee /etc/apt/sources.list.d/webupd8team-java.list \
	&& echo "deb-src http://ppa.launchpad.net/webupd8team/java/ubuntu xenial main" | tee -a /etc/apt/sources.list.d/webupd8team-java.list \
	&& apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys EEA14886 \
	&& apt-get update \
	&& apt-get install oracle-java8-installer

## Install igraph
RUN apt-get update && apt-get install -y --no-install-recommends \
	libigraph0 \
	libigraph0-dev

## Install other libraries to build rJava
RUN apt-get update && apt-get install -y --no-install-recommends \
	liblzma-dev \
	libbz2-dev \
	libicu-dev

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

RUN Rscript -e 'install.packages(pkgs = c("preText"), dependencies = TRUE,repos="https://mran.microsoft.com/snapshot/2016-12-25")'

# for john-henry
RUN Rscript -e 'install.packages(pkgs = c("rvest"), dependencies = TRUE,repos="https://mran.microsoft.com/snapshot/2016-12-25")'

#in-container data if folks don't want to volume map
RUN git clone https://github.com/dlabctawg/docks /home/rstudio

#CMD["\init"]
