# bioinformatics_scripts

## Preprocessing

To build a docker container do ./build_container.sh

To run docker container do ./run_container.sh

It emits a password for Rstudio, use it to log in on localhost:8787

When in Rstudio, run the script preproc.R

Once finished, use data_upload_script.ipynb to upload the data to s3

Run it in a folder with all the required files: `cells_output.tsv`, `genes_output.tsv`, `normalized.mtx`, `raw.mtx`.

## Data upload

Before you use the data_upload_script.ipynb, set up your aws config so it points to the right region:

    aws configure set default.region eu-west-1

Make sure that your aws credentials file has the correct credentials

