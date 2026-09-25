CONDA_ENV_DIR := $(shell dirname $(CONDA_EXE))
CONDA_ENV_NAME := beyond_secretome
ENV_FILE := environment.yml
LOCK_FILE := conda-lock.yml

.PHONY: all dependencies install-conda-lock lock dependencies-locked configure-env test-env clean

all:
	@echo "Available targets:"
	@echo "  make dependencies"
	@echo "  make lock"
	@echo "  make dependencies-locked"
	@echo "  make test-env"

dependencies:
	conda env update --name "$(CONDA_ENV_NAME)" --file "$(ENV_FILE)" --prune

install-conda-lock:
	"$(CONDA_EXE)" install --name base --channel conda-forge conda-lock

lock:
	"$(CONDA_ENV_DIR)/conda-lock" lock --file "$(ENV_FILE)" --platform linux-64 --platform osx-arm64 --platform osx-64 --lockfile "$(LOCK_FILE)"

dependencies-locked:
	"$(CONDA_ENV_DIR)/conda-lock" install --name "$(CONDA_ENV_NAME)" "$(LOCK_FILE)"
	$(MAKE) configure-env

configure-env:
	@env_prefix="$$( conda run --name "$(CONDA_ENV_NAME)" python -c 'import sys; print(sys.prefix)' )" && \
	echo "Setting R_LIBS_USER=$${env_prefix}/lib/R/library" && \
	conda env config vars set --name "$(CONDA_ENV_NAME)" R_LIBS_USER="$${env_prefix}/lib/R/library"

test-env:
	source "$(CONDA_ENV_DIR)/activate" "$(CONDA_ENV_NAME)" && \
	Rscript --vanilla -e 'packages <- c("readr", "readxl", "dplyr", "igraph", "ggplot2", "scales", "ggrepel", "Matrix", "ggraph", "tidygraph", "data.table", "tidyverse", "future", "furrr"); invisible(lapply(packages, library, character.only = TRUE)); cat("All R packages loaded successfully.\n")' && \
	python -c 'import gseapy as gp; from gseapy import enrichr; import pandas as pd; import seaborn as sns; import matplotlib.pyplot as plt; from sklearn.feature_extraction.text import TfidfVectorizer; import numpy as np; print("All Python packages loaded successfully.")'

clean:
	conda env remove --name "$(CONDA_ENV_NAME)"