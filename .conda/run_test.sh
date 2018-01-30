#!/bin/bash

set -ev

# Run tests inside the source directory!
cd "${RECIPE_DIR}/.."
flake8
python -m pytest --cov="${SP_DIR}/${PKG_NAME}"
