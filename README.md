# pH glider

Tools for processing and plotting glider data with a pH sensor on-board.

Author: Lori Garzio (lgarzio@marine.rutgers.edu)

**This repo is under development**

## Installation Instructions
Add the channel conda-forge to your .condarc. You can find out more about conda-forge from their website: https://conda-forge.org/

`conda config --add channels conda-forge`

Clone the phglider repository

`git clone https://github.com/lgarzio/phglider.git`

Change your current working directory to the location that you downloaded phglider. 

`cd /Users/garzio/Documents/repo/phglider/`

Create conda environment from the included environment.yml file:

`conda env create -f environment.yml`

Once the environment is done building, activate the environment:

`conda activate phglider`

Install the toolbox to the conda environment from the root directory of the phglider toolbox:

`pip install .`

The toolbox should now be installed to your conda environment.

## Glider Data Processing Steps
1. [download_ds.py](https://github.com/lgarzio/phglider/blob/master/delayed_analysis/download_ds.py): Download the delayed-mode dataset to your local machine.

2. Calculate the best time shifts for pH and dissolved oxygen using MATLAB code bestShifts_updown.m

3. [time_shift.py](https://github.com/lgarzio/phglider/blob/master/delayed_analysis/time_shift.py): Apply the best time shifts calculated in step 2.

4. [plot_time_shift.py](https://github.com/lgarzio/phglider/blob/master/delayed_analysis/plot_time_shift.py): Evaluate the time shifts by plotting short time ranges (e.g. 3 hour) for pH and dissolved oxygen.

5. [ta_sal_regression.py](https://github.com/lgarzio/phglider/blob/master/delayed_analysis/ta_sal_regression.py): Calculate the best fit TA-salinity linear regression when the discrete samples are analyzed.

6. [glider_qc.py](https://github.com/lgarzio/phglider/blob/master/delayed_analysis/glider_qc.py): Run QARTOD QA/QC tests.