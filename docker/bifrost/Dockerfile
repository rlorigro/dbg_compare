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

## Add time library for basic profiling
RUN apt-get install -y --no-install-recommends time

## Install bifrost from source
WORKDIR /software
RUN wget https://github.com/pmelsted/bifrost/archive/refs/tags/v1.2.0.tar.gz && tar -xvzf v1.2.0.tar.gz
WORKDIR /software/bifrost-1.2.0/build/
RUN cmake -DCOMPILATION_ARCH=x86-64 .. && make -j 4 && make install
ENV PATH=$PATH:/usr/local/bin/
ENV LIBRARY_PATH=$LIBRARY_PATH:/usr/local/lib/
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/

## Install python3
RUN apt-get install -y --no-install-recommends python3
RUN apt-get install -y --no-install-recommends python3-pip


## Download dbg_compare and force rebuild with time sensitive command
ADD http://date.jsontest.com /etc/builddate

ADD scripts/ /software

WORKDIR /software
