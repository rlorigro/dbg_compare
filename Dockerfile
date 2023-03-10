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

## Install samtools
WORKDIR /software
RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && tar -xvjf samtools-1.17.tar.bz2

WORKDIR samtools-1.17
RUN autoheader
RUN autoconf -Wno-syntax
RUN ./configure
RUN make -j 4

ENV PATH=/software/samtools-1.17:$PATH

## Install python3
RUN apt-get install -y --no-install-recommends python3
RUN apt-get install -y --no-install-recommends python3-pip

## Install some packages
RUN yes | pip3 install pysam
RUN yes | pip3 install google-auth
RUN yes | pip3 install requests

## Download dbg_compare and force rebuild with time sensitive command
ADD http://date.jsontest.com /etc/builddate

WORKDIR /software
RUN git clone https://github.com/rlorigro/dbg_compare.git

ENV PYTHONPATH="${PYTHONPATH}:/software/dbg_compare/scripts/"
