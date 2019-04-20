FROM meval/scns:latest

LABEL maintainer="Aditya Pratapa <adyprat@vt.edu>"

USER root

WORKDIR /

RUN pip2 install pandas==0.21

RUN git clone https://github.com/fionahamey/Pseudotime-network-inference

WORKDIR Pseudotime-network-inference/

RUN git checkout b900655

#RUN python2 booleanRules.py Bptf binary_expression_LMPP.txt 5 LMPP_trajectory_order.txt 2 2 startingNetworkParCor.txt 0.95 0.05

#RUN cat Bptf_boolean_rules_5.txt


RUN apt-get install time