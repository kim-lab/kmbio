#!/bin/bash

set -ex

if [[ \
      -z $1 ]] ; then
    echo "You must specify the test to run!"
    exit -1
fi

mkdir -p conda-bld public

CI_PROJECT_NAME="$(basename $(pwd))"
CI_PROJECT_NAMESPACE="$(basename $(dirname $(pwd)))"

echo "${CI_PROJECT_NAME}"

gitlab-runner exec docker \
    --env CI_PROJECT_NAME="${CI_PROJECT_NAME}" \
    --env HOST_USER_ID="$(id -u $USER)" \
    --docker-volumes "$(pwd):/home/${CI_PROJECT_NAME}:rw" \
    $1
