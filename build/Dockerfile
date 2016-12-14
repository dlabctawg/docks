FROM rocker/tidyverse:3.3.2

RUN Rscript -e 'install.packages(pkgs = c("igraph","data.table","formatR","stm","stmBrowser","stmCorrViz","LDAvis","stringdist","qdap","qdapDictionaries","mallet","RTextTools"), dependencies = TRUE)'

RUN apt-get update && apt-get install -y --no-install-recommends \
	default-jdk

RUN R CMD javareconf

RUN Rscript -e 'devtools::install_github("brooksambrose/JSTORr")'

RUN cd ~ && git clone https://github.com/dlabctawg/jstorr-demo && unzip jstorr-demo/*zip

#CMD["\init"]
