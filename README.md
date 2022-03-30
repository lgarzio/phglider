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

## Delayed Mode Glider Data Processing Steps
1. [download_ds.py](https://github.com/lgarzio/phglider/blob/master/delayed_analysis/download_ds.py): Download the delayed-mode dataset to your local machine from the [RUCOOL ERDDAP](http://slocum-data.marine.rutgers.edu/erddap/index.html) server.

2. [ta_sal_regression.py](https://github.com/lgarzio/phglider/blob/master/ta_equation/ta_sal_regression.py): Calculate the best fit TA-salinity linear regression when the discrete samples are analyzed.

3. [process_phglider.py](https://github.com/lgarzio/phglider/blob/master/delayed_analysis/process_phglider.py): Process delayed-mode pH glider data.
   1. Apply QARTOD QC flags to CTD data.
   2. Apply CTD hysteresis test flags.
   3. Run location QARTOD test.
   4. Convert pH voltages of 0.0 to nan.
   5. Interpolate (method=linear) pressure and remove data where pressure <1 dbar.
   6. Interpolate (method=linear) temperature and salinity.
   7. Calculate pH from original and corrected voltages and interpolated CTD data.
   8. Calculate Total Alkalinity from interpolated salinity using a linear relationship determined from in-situ water sampling.
   9. Run ioos_qc gross range test on additional variables defined in the gross_range.yml config file, and apply test results to data.
   10. Run QARTOD spike test on pH and corrected pH, and apply test results to data.
   11. Calculate CO2SYS variables using TA, corrected pH, interpolated salinity, interpolated temperature, interpolated pressure.
   12. Convert oxygen concentration to mg/L.

4. [plot_qc.py](https://github.com/lgarzio/phglider/blob/master/delayed_analysis/plot_qc.py): Create cross-section and profile plots for short sections of glider data for the purposes of QC.

5. Create [configuration files](https://github.com/lgarzio/phglider/tree/master/config) for the deployment.

6. [glider_to_dac.py](https://github.com/lgarzio/phglider/blob/master/delayed_analysis/glider_to_dac.py): Format dataset to upload to the [IOOS glider DAC](https://gliders.ioos.us/).

## Citations
Humphreys, M. P., Gregor, L., Pierrot, D., van Heuven, S. M. A. C., Lewis, E. R., and Wallace, D. W. R. (2020). [PyCO2SYS](https://pypi.org/project/PyCO2SYS/): marine carbonate system calculations in Python. Zenodo. doi:10.5281/zenodo.3744275.

Lewis, E. and Wallace, D. W. R. (1998) Program Developed for CO2 System Calculations, ORNL/CDIAC-105, Carbon Dioxide Inf. Anal. Cent., Oak Ridge Natl. Lab., Oak Ridge, Tenn., 38 pp., [https://salish-sea.pnnl.gov/media/ORNL-CDIAC-105.pdf](https://salish-sea.pnnl.gov/media/ORNL-CDIAC-105.pdf)

Saba, G.K., Wright‐Fairbanks, E., Chen, B., Cai, W.‐J., Barnard, A. H., Jones, C. P., et al. (2019). The development and validation of a profiling glider deep ISFET‐based pH sensor for high resolution observations of coastal and ocean acidification. Frontiers in Marine Science, 6. https://doi.org/10.3389/fmars.2019.00664