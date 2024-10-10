<p align="center">
  <img src="image/logo.png"/>
</p>

---

Atom-based machine learning for estimating nucleophilicity and electrophilicity by predicting methyl cation affinities (MCAs) and methyl anion affinities
(MAAs), respectively.

[TRY IT HERE: https://www.esnuel.org](https://www.esnuel.org)

GitHub repository for our QM-based workflow is: https://github.com/jensengroup/esnuel

GitHub repository for our atom-based descriptor: https://github.com/NicolaiRee/smi2gcs


## Installation

Clone this Github repository:

    git clone https://github.com/jensengroup/ESNUEL_ML.git

Download and unpack models:

    wget -O models.tar.xz https://sid.erda.dk/share_redirect/Ear6g6wl0G; tar xf models.tar.xz; mv models src/esnuelML/.

Download and unpack datasets (OBS! Not required for running the models):

    wget -O data.tar.xz https://sid.erda.dk/share_redirect/c7LF5NaYvH; tar xf data.tar.xz

Install the Python environment using `conda`:

    conda env create -f etc/environment.yml && conda activate esnuelML

Finally, download the binaries of xtb version 6.7.0:

    mkdir dep; cd dep; wget https://github.com/grimme-lab/xtb/releases/download/v6.7.0/xtb-6.7.0-linux-x86_64.tar.xz; tar -xvf ./xtb-6.7.0-linux-x86_64.tar.xz; cd ..


## Usage

An example of usage via CLI command:

    python src/esnuelML/predictor.py --smi 'CCOC(=O)CCN(SN(C)C(=O)Oc1cccc2c1OC(C)(C)C2)C(C)C'

The results are saved in a subfolder under "./desc_calcs" with a visual output of the predictions (saved in a .html file).


## Raw QM calculation results

All the raw input and output files from the QM calculations can be downloaded here (compressed size: 29 GB, unpacked size: 355 GB):

    wget -O calculations.tar.xz https://sid.erda.dk/share_redirect/hp3ttQcTgs; tar xf calculations.tar.xz


## Citation 

Our work is available as a preprint on [ChemRxiv](https://doi.org/), where more information is available.
```
@article{ree2024esnuelML,
  title = {Atom-Based Machine Learning for Estimating Nucleophilicity and Electrophilicity with Applications to Retrosynthesis and Chemical Stability},
  url = {http://dx.doi.org/},
  DOI = {},
  author = {Nicolai Ree and Jan M. Wollschl\"{a}ger and Andreas H. G\"{o}ller and Jan H. Jensen},
  year = {2024}
}
```
