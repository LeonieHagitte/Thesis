FROM rocker/verse:4.2.3
WORKDIR /home/rstudio
RUN apt-get update -y && apt-get install -y rsync
RUN tlmgr install collection-latexrecommended
