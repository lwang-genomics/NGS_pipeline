#!/bin/bash
set -e

mkdir lib
cd lib

echo "[INFO] Downloading genome and index files..."

URL1="https://download_link.com/lib_hg38.tar.gz"
URL2="https://download_link.com/lib_mm10.tar.gz"

# Download and extract
wget -O hg38.tar.gz $URL1
tar -xzf hg38.tar.gz

wget -O mm10.tar.gz $URL2
tar -xzf mm10.tar.gz

rm hg38.tar.gz mm10.tar.gz

echo "[INFO] Genome references extracted into ./lib/"
