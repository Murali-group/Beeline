# ModelEval

This is the main repository for GRN model evaluation

Instructions
- To install docker on Ubuntu 18.04, follow the steps mentioned [here](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-18-04)
- Setup docker to run docker without sudo using ` sudo usermod -aG docker $USER`, if you haven't already. See more details [here](https://askubuntu.com/questions/477551/how-can-i-use-docker-without-sudo).
- To create the docker containers for each of the algorithms run `. initialize.sh`
- To compute ranked list of edges, run `python eval.py --config config-files/config.yaml`
- To run the evaluation pipeline, run `python eval_summarizer.py --config config-files/config.yaml`
