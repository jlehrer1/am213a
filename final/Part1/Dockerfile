FROM ubuntu

WORKDIR /src

RUN apt-get update
RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN apt-get --allow-releaseinfo-change update && \
    apt-get install -y --no-install-recommends \
        curl \
        sudo \
        vim \ 
        build-essential \ 
        libblas-dev \
        liblapack-dev \
        gfortran

COPY . .