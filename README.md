# ModelEval

This is the main repository for GRN model evaluation

Instructions

- Setup docker to run docker without sudo using ` sudo usermod -aG docker $USER`, if you haven't already. See more details [here](https://askubuntu.com/questions/477551/how-can-i-use-docker-without-sudo).
- To create the docker containers for each of the algorithms run `. initialize.sh`
- To run te evaluation pipeline, run `python eval.py --config config/config.yaml`
