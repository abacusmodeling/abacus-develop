FROM debian:bullseye-slim

RUN apt-get update && apt-get install -y --no-install-recommends git gfortran libssl-dev make cmake vim wget bc \
    && apt-get install -y --no-install-recommends mpich libmpich-dev

ENV GIT_SSL_NO_VERIFY 1

RUN export cmake_url="https://github.com/Kitware/CMake/releases/download/v3.20.0/cmake-3.20.0-linux-x86_64.sh" \
    && wget --no-check-certificate --quiet $cmake_url \
    && export file=$(basename "$cmake_url") \
    && bash $file --prefix=/usr --skip-license \
    && rm $file

RUN cd /tmp \
    && git clone https://github.com/USCiLab/cereal.git \
    && cp -r cereal/include /usr/local \
    && rm -rf cereal

RUN cd /tmp \
    && git clone https://github.com/xianyi/OpenBLAS.git --single-branch --depth=1 \
    && cd OpenBLAS && mkdir build && cd build \
    && cmake .. \
    && cmake --build . -j9 \
    && cmake --install . --prefix /opt \
    && cd /tmp && rm -rf OpenBLAS

RUN cd /tmp \
    && git clone https://github.com/darelbeida/scalapack.git -b v2.0.2-openblas --single-branch --depth=1 \
    && cd scalapack && make lib && cp libscalapack.a /usr/local/lib/ \
    && cd /tmp && rm -rf scalapack

RUN cd /tmp \
    && git clone https://github.com/darelbeida/elpa.git -b ELPA_2016.05.004-openblas --single-branch --depth=1 \
    && cd elpa && mkdir build && cd build \
    && ../configure CFLAGS="-O3 -march=native -mavx2 -mfma -funsafe-loop-optimizations -funsafe-math-optimizations -ftree-vect-loop-version -ftree-vectorize" \
        FCFLAGS="-O2 -mavx" \
    && make -j8 && make PREFIX=/usr/local install \
    && ln -s /usr/local/include/elpa-2016.05.004/elpa /usr/local/include/ \
    && cd /tmp && rm -rf elpa

RUN cd /tmp \
    && wget http://www.fftw.org/fftw-3.3.9.tar.gz \
    && tar zxvf fftw-3.3.9.tar.gz \
    && cd fftw-3.3.9 \
    && ./configure --enable-mpi-fortran --enable-orterun-prefix-by-default FC=gfortran \
    && make -j8 && make PREFIX=/usr/local install \
    && cd /tmp && rm -rf fftw-3.3.9 && rm fftw-3.3.9.tar.gz

ENV LD_LIBRARY_PATH /usr/local/lib

RUN apt-get install -y unzip

RUN cd /tmp \
    && wget https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-1.9.0%2Bcpu.zip --no-check-certificate --quiet \
    && unzip libtorch-shared-with-deps-1.9.0+cpu.zip \
    && cp -r libtorch/include /usr/local \
    && cp -r libtorch/lib /usr/local \
    && cp -r libtorch/share /usr/local \
    && rm -rf libtorch

RUN cd /tmp \
    && wget https://gitlab.com/libxc/libxc/-/archive/5.1.5/libxc-5.1.5.tar.gz --no-check-certificate --quiet \
    && tar xzf libxc-5.1.5.tar.gz \
    && cd libxc-5.1.5 \
    && mkdir build \
    && cmake -B build -DBUILD_TESTING=OFF \
    && cmake --build build \
    && cmake --install build \
    && cd /tmp \
    && rm -rf libxc-5.1.5 \
    && rm libxc-5.1.5.tar.gz
