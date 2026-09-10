#!/usr/bin/env bash
# Idempotent development-environment bootstrap for the CFPerm repository.
# Installs the R toolchain (R + grf + testthat + devtools), the Python analysis
# dependencies, and unpacks the bundled datasets used by the plotting scripts.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "== [1/4] System packages (R, build toolchain, dev libraries) =="
export DEBIAN_FRONTEND=noninteractive
sudo apt-get update -qq
sudo apt-get install -y -qq --no-install-recommends \
  r-base r-base-dev build-essential gfortran \
  libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
  r-cran-testthat r-cran-devtools r-cran-roxygen2 r-cran-mass \
  python3 python3-pip unzip

echo "== [2/4] R package: grf (compiled from CRAN if missing) =="
if ! Rscript -e 'quit(status = !requireNamespace("grf", quietly = TRUE))'; then
  sudo Rscript -e 'install.packages("grf", repos = "https://cloud.r-project.org", Ncpus = parallel::detectCores())'
fi
Rscript -e 'cat("grf", as.character(packageVersion("grf")), "\n")'

echo "== [3/4] Python dependencies =="
python3 -m pip install --break-system-packages -q -r requirements.txt

echo "== [4/4] Unpack bundled datasets =="
mkdir -p datasets/extracted
unzip -o -j datasets/datasets.zip \
  'source_DiabetesReadmission.csv' 'target_DiabetesReadmission.csv' \
  -d datasets/extracted >/dev/null

echo "Environment setup complete."
