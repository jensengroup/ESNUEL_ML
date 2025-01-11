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


## Model training
The scripts (.py files) to train/retrain the ML models are located in the following directory: scripts/MLtrain/

    - elec_SMI2GCS_3_cm5: The lightGBM model used to predict the MAA values, referred to as the MAA ML model in the paper.
    - elec_RF_SMI2GCS_3_cm5: The random forest model used for uncertainty quantification of the MAA ML model predictions.
    - nuc_SMI2GCS_3_cm5: The lightGBM model used to predict the MCA values, referred to as the MCA ML model in the paper.
    - nuc_RF_SMI2GCS_3_cm5: The random forest model used for uncertainty quantification of the MCA ML model predictions.

Before running the training scripts (.py files) in the above folders, please make sure to download and unpack datasets:

    wget -O data.tar.xz https://sid.erda.dk/share_redirect/c7LF5NaYvH; tar xf data.tar.xz

Then change the path to the data in all training scripts just below the headline "Training/test data split":
    
    df = pd.read_csv('insert_full_path_to_esnuelML/data/df_elec_x.csv.gz', index_col=0)

You are now ready to run the training scripts e.g.:

    python LGBMRegressor_optuna.py

Once the training is completed the lightGBM models are saved as "final_best_model.txt" and the radom forest models are saved as "final_best_model.joblib" in the respective folders.
We then recommend to compress the lightGBM models using the provided bash script "do_gzip.sh" by simply running the script:

    bash do_gzip.sh

Hereafter, rename the model files and copy them into the correct model folders: 

    scripts/MLtrain/elec_SMI2GCS_3_cm5/SMI2GCS_3_cm5_model.txt.gz -> src/esnuelML/models/elec/SMI2GCS_3_cm5_model.txt.gz
    scripts/MLtrain/elec_RF_SMI2GCS_3_cm5/final_best_model.joblib -> src/esnuelML/models/elec/RF_SMI2GCS_3_cm5.joblib
    scripts/MLtrain/nuc_SMI2GCS_3_cm5/SMI2GCS_3_cm5_model.txt.gz -> src/esnuelML/models/nuc/SMI2GCS_3_cm5_model.txt.gz
    scripts/MLtrain/nuc_RF_SMI2GCS_3_cm5/final_best_model.joblib -> src/esnuelML/models/nuc/RF_SMI2GCS_3_cm5.joblib


## Citation 

Our work is available as a preprint on [ChemRxiv](https://doi.org/10.26434/chemrxiv-2024-2p7ch), where more information is available.
```
@article{ree2024esnuelML,
  title = {Atom-Based Machine Learning for Estimating Nucleophilicity and Electrophilicity with Applications to Retrosynthesis and Chemical Stability},
  url = {http://dx.doi.org/10.26434/chemrxiv-2024-2p7ch},
  DOI = {10.26434/chemrxiv-2024-2p7ch},
  author = {Ree,  Nicolai and Wollschl\"{a}ger,  Jan M. and G\"{o}ller,  Andreas H. and Jensen,  Jan H.},
  year = {2024},
  month = oct 
}
```
