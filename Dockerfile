# syntax = docker/dockerfile:experimental


FROM python:3.12.3


RUN --mount=type=cache,target=/var/cache/apt \
    echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections \
    && apt-get update -y \
    && apt-get install -y --no-install-recommends \
        build-essential \
        pkg-config \
        python3-dev \
        libhdf5-dev

# python 3.12 may have issues w/ pip:
# AttributeError: module 'pkgutil' has no attribute 'ImpImporter'. Did you mean: 'zipimporter'?
RUN python -m ensurepip --upgrade \
    && python -m pip install --upgrade setuptools

WORKDIR /tmp
COPY setup.py .
COPY process_genome.py .
RUN mkdir transposon
COPY transposon/* transposon/
RUN --mount=type=cache,mode=0755,target=/root/.cache/pip \
    python -m pip install .
