FROM pkharchenkolab/pagoda2:latest

RUN mkdir -p /home/rstudio/script

COPY preproc.R /home/rstudio/script
