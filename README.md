# bioinformatics_scripts

preproc.R is best run in a docker container:

    docker run -e PASSWORD=<your password> -p 8787:8787 -v /Users/adamkurkiewicz/path_to_files:/home/rstudio/path_to_files pkharchenkolab/pagoda2:latest

# upload.py

Set up your aws config so it points to the right region:

    aws configure set default.region eu-west-1

Then you can run `upload.py` in a folder with all the required files: `cells_output.tsv`, `genes_output.tsv`, `normalized.mtx`, `raw.mtx`.