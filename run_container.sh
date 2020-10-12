password=$(base64 < /dev/urandom | tr -d 'O0Il1+/' | head -c 8)

echo "Log in to R session on localhost:8787 with password=$password"

path_to_experiment="/Users/adamkurkiewicz/programming/ange_control"

docker run -e PASSWORD="$password" -p 8787:8787 -v "$path_to_experiment:/home/rstudio/data" biomage_scriptrunner

