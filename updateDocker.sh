## This docker command updates the blast packages required. It builds over
## the image original VirSorter docker image instead of rebuilding from 
## scratch

docker build -f ./Dockerfile-updateBLAST -t discoenv/virsorter:updateBLAST .
