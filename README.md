# bioinformatics_scripts

preproc.R is best run in a docker container:


docker run -e PASSWORD=<your password> -p 8787:8787 -v /Users/adamkurkiewicz/path_to_files:/home/rstudio/path_to_files pkharchenkolab/pagoda2:latest
