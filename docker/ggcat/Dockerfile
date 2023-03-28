FROM ubuntu:22.04

MAINTAINER rlorigro@broadinstitute.edu

ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

## Foundational
RUN apt-get update
RUN apt-get install -y --reinstall ca-certificates
RUN update-ca-certificates -f

## Basic software
RUN apt-get install -y --no-install-recommends wget
RUN apt-get install -y --no-install-recommends curl
RUN apt-get install -y --no-install-recommends cmake
RUN apt-get install -y --no-install-recommends autoconf
RUN apt-get install -y --no-install-recommends build-essential
RUN apt-get install -y --no-install-recommends bzip2
RUN apt-get install -y --no-install-recommends git
RUN apt-get install -y --no-install-recommends sudo
RUN apt-get install -y --no-install-recommends pkg-config
RUN apt-get install -y --no-install-recommends zlib1g-dev
RUN apt-get install -y --no-install-recommends libcurl4-openssl-dev
RUN apt-get install -y --no-install-recommends libbz2-dev
RUN apt-get install -y --no-install-recommends libncurses5-dev
RUN apt-get install -y --no-install-recommends liblzma-dev

## Install rust
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH=$PATH:/root/.cargo/bin
RUN echo $PATH
RUN rustup update

## Install ggcat from source
WORKDIR /software
RUN git clone https://github.com/algbio/ggcat --recursive && cd ggcat && git checkout a256dd5
WORKDIR /software/ggcat/
RUN cargo install --path crates/cmdline/ --locked

## Add time library for basic profiling
RUN apt-get install -y --no-install-recommends time

## Install python3
RUN apt-get install -y --no-install-recommends python3
RUN apt-get install -y --no-install-recommends python3-pip

## Download dbg_compare and force rebuild with time sensitive command
ADD http://date.jsontest.com /etc/builddate

ADD scripts/ /software

WORKDIR /software
