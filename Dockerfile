FROM continuumio/miniconda3:4.7.12
ARG fname=AmberTools19.tar.bz2
SHELL ["/bin/bash", "-c"]
ENV AMBERHOME=/rism/amber18
ENV LOG_LEVEL=info

WORKDIR /rism

RUN wget "http://ambermd.org/cgi-bin/AmberTools19-get.pl?Name=bot&Institution=NA&City=NA&State=NA&Country=NA&OS=linux64" -O $fname
#COPY ../$fname $fname

RUN tar -xvf $fname && \
  apt-get update --fix-missing && \
  apt-get install -y bc csh flex gfortran g++ xorg-dev zlib1g-dev libbz2-dev patch bison build-essential && \
  cd $AMBERHOME && \
  echo y | ./configure -noX11 --with-python $(which python) --python-install global gnu && \
  source $AMBERHOME/amber.sh && make install

RUN echo "source /rism/amber18/amber.sh" >> ~/.bashrc

WORKDIR /python-rism
COPY ./ ./

RUN chmod +x run-gunicorn.sh
RUN conda env create -f environment.yml

# Make RUN commands use the new environment:
RUN echo "source activate syntelly-calc" >> ~/.bashrc

ENV PATH /opt/conda/envs/syntelly-calc/bin:$PATH

SHELL ["conda", "run", "-n", "syntelly-calc", "/bin/bash", "-c"]

RUN pip install gunicorn[gevent]

EXPOSE 3002

ENTRYPOINT ["./run-gunicorn.sh"]
