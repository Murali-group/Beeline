FROM meval/modeleval:latest

LABEL maintainer="Aditya Pratapa <adyprat@vt.edu>"

USER root
ENV TZ=America/New_York
ENV DEBIAN_FRONTEND=noninteractive 
RUN echo "deb http://archive.debian.org/debian stretch main contrib non-free" > /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 04EE7237B7D453EC
#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 648ACFD622F3D138
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 0E98404D386FA1D9
#RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys EF0F382A1A7B6500
RUN apt-get update --allow-unauthenticated
RUN apt-get install -y build-essential patch --allow-unauthenticated
RUN apt-get install -y fsharp \ 
	mono-xbuild --allow-unauthenticated
RUN git clone https://github.com/Z3Prover/z3
WORKDIR z3
RUN git checkout d6df51951f4cdc95f0dfd3b1297d04a465d8f2ca 

RUN python2 scripts/mk_make.py --python
RUN cd build && make
RUN cd build && make install
WORKDIR ..

RUN export PATH="$PATH:/z3/build/"

RUN  git clone https://github.com/swoodhouse/SCNS-Toolkit
WORKDIR SCNS-Toolkit
RUN git checkout 27cb7a349f239d450a9571b270abc38b053ad6b2

WORKDIR SynthesisEngine

RUN xbuild SynthesisEngine.sln

RUN cp bin/Release/* .

#RUN mkdir outDir
#RUN mono SynthesisEngine.exe cmpStates.csv cmpEdges.csv cmpParameters.csv cmp_initial_states.txt cmp_target_states.txt outDir
#WORKDIR outDir/
#RUN cat Cebpa.txt

RUN apt-get install time


