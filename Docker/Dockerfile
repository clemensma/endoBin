FROM ubuntu:22.04

LABEL author="Clemens Mauksch"

RUN apt update --fix-missing && apt install sudo
RUN sudo apt update --fix-missing

RUN sudo apt install -y apt-utils python3 cutadapt fastqc wget git ncbi-blast+ trim-galore make gcc zlib1g-dev --fix-missing
RUN sudo apt install -y pigz

RUN adduser --disabled-password --gecos '' endominer
RUN adduser endominer sudo
RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
USER endominer
WORKDIR /home/endominer/Docker
RUN sudo  chmod a+rwx /home/endominer/Docker

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b
RUN rm Miniconda3-latest-Linux-x86_64.sh 

ENV PATH /home/endominer/miniconda3/bin:$PATH

RUN conda update conda
RUN conda update --all

RUN conda install -c anaconda pip
RUN pip install better_fasta_grep
RUN conda install -c bioconda seqkit

RUN git clone https://github.com/lh3/seqtk.git;
RUN cd seqtk; make