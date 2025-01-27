# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# ## boilerplate

FROM ubuntu:24.04

ENV TZ=UTC
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get upgrade -y && \
    apt-get install -y --no-install-recommends --fix-missing \
    apt-utils \
    build-essential \
    bash \
    wget \
    git \
    curl \
    tar \
    unzip \
    patch \
    gzip \
    bzip2 \
    file \
    gfortran \
    coreutils \
    libopenblas-dev \
    liblapack-dev \
    pkg-config \
    automake \
    make \
    cmake \
    gnupg \
    cmake-curses-gui \
    libreadline-dev \
    libxml2-dev \
    software-properties-common \
    python3.10 python3-dev python3-pip \
    python3-setuptools python3-numpy python3-scipy \
    libeigen3-dev libyaml-cpp-dev libopenmpi-dev libparmetis-dev libmetis-dev ninja-build && \
    wget -O llvm-snapshot.gpg.key https://apt.llvm.org/llvm-snapshot.gpg.key && \
    apt-key add llvm-snapshot.gpg.key && \
    rm -f llvm-snapshot.gpg.key && \
    add-apt-repository "deb http://apt.llvm.org/noble/ llvm-toolchain-noble-19 main" && \
    add-apt-repository "deb-src http://apt.llvm.org/noble/ llvm-toolchain-noble-19 main" && \
    apt-get -y update && \
    apt-get -y install clang-19 clang-tidy-19 libomp-19-dev

# ## IO

# HDF5
RUN wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.3/src/hdf5-1.12.3.tar.bz2 && \
    tar -xf hdf5-1.12.3.tar.bz2 && cd hdf5-1.12.3 && \
    CPPFLAGS="-fPIC ${CPPFLAGS}" CC=mpicc CXX=mpicxx ./configure --enable-parallel --prefix=/usr --with-zlib --disable-shared && \
    make -j $(nproc) install && \
    cd .. && rm -rf hdf5-1.12.3

# Netcdf
RUN wget https://downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz && \
    tar -xf netcdf-c-4.9.2.tar.gz && cd netcdf-c-4.9.2 && \
    CFLAGS="-fPIC ${CFLAGS}" CC=h5pcc ./configure --enable-shared=no --prefix=/usr --disable-dap --disable-byterange && \
    make -j $(nproc) install && \
    cd .. && rm -rf netcdf-c-4.9.2

# ## input data

# ImpalaJIT
RUN git clone --depth 1 https://github.com/uphoffc/ImpalaJIT.git gitbuild && \
    mkdir -p gitbuild/build && cd gitbuild/build && \
    cmake .. -GNinja && \
    ninja install && \
    cd ../.. && rm -rf gitbuild

# ASAGI
RUN git clone --recursive --depth 1 https://github.com/TUM-I5/ASAGI.git gitbuild && \
mkdir -p gitbuild/build && cd gitbuild/build && \
    cmake .. -DSHARED_LIB=OFF -DFORTRAN_SUPPORT=OFF -DSTATIC_LIB=ON -GNinja && \
    ninja install && \
    cd ../.. && rm -rf gitbuild

# Lua
RUN wget https://www.lua.org/ftp/lua-5.3.6.tar.gz && \
    tar -xzvf lua-5.3.6.tar.gz && cd lua-5.3.6/ &&\
    make linux && make local && \
    cp -r install/* /usr && \
    cd .. && rm -rf lua-5.3.6

# easi
RUN git clone --recursive --depth 1 --branch v1.4.0 https://github.com/SeisSol/easi.git gitbuild && \
    mkdir -p gitbuild/build && cd gitbuild/build && \
    cmake .. -GNinja -DLUA=ON -DIMPALAJIT=ON -DASAGI=ON && \
    ninja install && \
    cd ../.. && rm -rf gitbuild

# ## code generators

# libxsmm
RUN git clone --depth 1 --branch 1.17 https://github.com/libxsmm/libxsmm.git gitbuild && \
    mkdir -p gitbuild && cd gitbuild && \
    make generator -j $(nproc) && cp bin/libxsmm_gemm_generator /usr/bin && \
    cd .. && rm -rf gitbuild

# PSpaMM
RUN pip3 install --break-system-packages git+https://github.com/SeisSol/PSpaMM.git

# ## setup and cleanup

WORKDIR /app

COPY . /app

RUN chmod +x /usr/bin/*

RUN apt-get clean && rm -rf /var/lib/apt/lists/*
