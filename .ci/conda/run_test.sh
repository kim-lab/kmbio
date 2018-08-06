#!/bin/bash

set -ev

python -m pytest \
    -c setup.cfg \
    --cov="${SP_DIR}/kmbio" \
    --benchmark-disable \
    --color=yes \
    tests/
