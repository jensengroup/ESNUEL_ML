<p align="center">
  <img src="image/logo.png"/>
</p>

---

Estimating Nucleophilicity and Electrophilicity From Atom-Based Machine Learning Predictions of Methyl Affinities


## Installation

For the installation, we recommend using `conda` to get all the necessary dependencies:

    conda env create -f environment.yml && conda activate esnuelML


Then download the binaries of xtb version 6.7.0:

    mkdir dep; cd dep; wget https://github.com/grimme-lab/xtb/releases/download/v6.7.0/xtb-6.7.0-linux-x86_64.tar.xz; tar -xvf ./xtb-6.7.0-linux-x86_64.tar.xz; cd ..


## Usage

An example of usage via CLI command:

    python run_ester_predictions.py --smi 'CCOC(=O)CCN(SN(C)C(=O)Oc1cccc2c1OC(C)(C)C2)C(C)C'

The xtb calculations are now saved in a "./calculations" folder.