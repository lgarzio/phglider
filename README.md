# pH glider

Tools for processing and plotting glider data with a pH sensor on-board.

Author: Lori Garzio (lgarzio@marine.rutgers.edu)

**_NOTE: this repo is under development_**

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

5. [ta_sal_regression.py](https://github.com/lgarzio/phglider/blob/master/ta_equation/ta_sal_regression.py): Calculate the best fit TA-salinity linear regression when the discrete samples are analyzed.

6. [glider_proc_qc.py](https://github.com/lgarzio/phglider/blob/master/delayed_analysis/glider_proc_qc.py): Calculate TA and CO2SYS variables, and run QARTOD QA/QC tests.

### Citations
Humphreys, M. P., Gregor, L., Pierrot, D., van Heuven, S. M. A. C., Lewis, E. R., and Wallace, D. W. R. (2020). [PyCO2SYS](https://pypi.org/project/PyCO2SYS/): marine carbonate system calculations in Python. Zenodo. doi:10.5281/zenodo.3744275.

Lewis, E. and Wallace, D. W. R. (1998) Program Developed for CO2 System Calculations, ORNL/CDIAC-105, Carbon Dioxide Inf. Anal. Cent., Oak Ridge Natl. Lab., Oak Ridge, Tenn., 38 pp., [https://salish-sea.pnnl.gov/media/ORNL-CDIAC-105.pdf](https://salish-sea.pnnl.gov/media/ORNL-CDIAC-105.pdf)