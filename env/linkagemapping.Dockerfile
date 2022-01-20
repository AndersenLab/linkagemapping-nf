FROM continuumio/miniconda
MAINTAINER Katie Evans <kathryn.evans@northwestern.edu>

RUN conda install R=3.6.0
RUN Rscript -e "install.packages('mscorefonts', dependencies = TRUE, repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('tidyverse', dependencies = TRUE, repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('devtools', dependencies = TRUE, repos = 'http://cran.us.r-project.org')"
RUN Rscript -e "install.packages('qtl', dependencies = TRUE, repos = 'http://cran.us.r-project.org')"
RUN conda install -c conda-forge tar
RUN Rscript -e "devtools::install_github('andersenlab/linkagemapping', dependencies = TRUE, repos = 'http://cran.us.r-project.org')"
RUN apt-get --allow-releaseinfo-change update && apt-get install -y procps  
