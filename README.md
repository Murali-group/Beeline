# :honeybee: BEELINE: Benchmarking gEnE reguLatory network Inference from siNgle-cEll transcriptomic data :honeybee:

This is the main repository for BEELINE. The documentation is available at: [https://murali-group.github.io/Beeline/](https://murali-group.github.io/Beeline/).

Instructions
- To install docker on Ubuntu 18.04, follow the steps mentioned [here](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-18-04)
- Setup docker to run docker without sudo using ` sudo usermod -aG docker $USER`, if you haven't already. See more details [here](https://askubuntu.com/questions/477551/how-can-i-use-docker-without-sudo).
- To create the docker containers for each of the algorithms run `sudo bash initialize.sh` (this step will take a while)
- To compute ranked list of edges, run `python eval.py --config config-files/config.yaml`
- To compute areas under the ROC and PR curves using the BEELINE's evaluation pipeline, run `python BLEvalAggregator.py --config config-files/config.yaml --auc`. To display the complete list of evalutation options, run `python BLEvalAggregator.py --help`.


If you use BEELINE in your research, please cite:

Aditya Pratapa, Amogh Jalihal, Jeffrey Law, Aditya Bharadwaj, and T M Murali. [Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data](https://doi.org/10.1101/642926), _bioRxiv_ (2019). doi.org/10.1101/642926

The repository for BoolODE is located at: [https://github.com/Murali-group/BoolODE](https://github.com/Murali-group/BoolODE)
