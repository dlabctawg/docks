FROM rocker/tidyverse:3.3.2

RUN Rscript -e 'install.packages(pkgs = c("igraph","data.table"), dependencies = TRUE)'

RUN apt-get update && apt-get install -y --no-install-recommends \
	default-jdk

RUN R CMD javareconf

RUN Rscript -e 'devtools::install_github("benmarwick/JSTORr")'

RUN cd ~ && git clone https://github.com/dlabctawg/jstorr-demo

RUN cd ~ && unzip jstorr-demo/*.zip

RUN Rscript -e 'install.packages(pkgs = c("stm"), dependencies = TRUE)'

CMD["/init"]
