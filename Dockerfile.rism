FROM continuumio/miniconda3:4.7.12
ARG fname=AmberTools19.tar.bz2
SHELL ["/bin/bash", "-c"]
ENV AMBERHOME=/source/amber18


WORKDIR /source

#RUN wget "http://ambermd.org/cgi-bin/AmberTools19-get.pl?Name=bot&Institution=NA&City=NA&State=NA&Country=NA&OS=linux64" -O $fname
COPY ./$fname $fname

RUN tar -xvf $fname && \
  apt-get update --fix-missing && \
  apt-get install -y bc csh flex gfortran g++ xorg-dev zlib1g-dev libbz2-dev patch bison build-essential && \
  cd $AMBERHOME && \
  echo y | ./configure -noX11 --with-python $(which python) --python-install global gnu && \
  source $AMBERHOME/amber.sh && make install


WORKDIR /opt

RUN conda install --yes numpy tini && \
  git clone https://github.com/MTS-Strathclyde/PC_plus pcplus


ENTRYPOINT ["tini", "-g", "--"]

