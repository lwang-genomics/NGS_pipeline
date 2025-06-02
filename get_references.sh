#!/bin/bash
set -e

mkdir -p  lib
cd lib

echo "[INFO] Downloading genome and index files..."

URL1="https://zenodo.org/records/15571133/files/hg38.tar.gz?download=1"
URL2="https://zenodo.org/records/15571143/files/mm10.tar.gz?download=1"

# Download and extract
wget -O hg38.tar.gz $URL1
tar -xzf hg38.tar.gz

wget -O mm10.tar.gz $URL2
tar -xzf mm10.tar.gz

rm hg38.tar.gz mm10.tar.gz

cd ..

echo "[INFO] Genome references extracted into ./lib/"
