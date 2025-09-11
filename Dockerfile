#FROM intel/fortran-essentials
#FROM intel/oneapi-hpckit:2022.1.1-devel-ubuntu18.04

FROM intel/oneapi-hpckit:latest

# Set environment variables for noninteractive install
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
 RUN apt-get update && apt-get install -y \
    wget \
    gnupg \
    make \
    build-essential \
    bash \
    cmake \ 
    git \
    && rm -rf /var/lib/apt/lists/*

# Always source Intel environment
 RUN echo "source /opt/intel/oneapi/setvars.sh" >> /etc/bash.bashrc --force

# Default shell
SHELL ["/bin/bash", "-c"]

CMD ["bash"] 
