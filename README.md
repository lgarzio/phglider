# pH glider

Tools for processing and plotting glider data with a pH sensor on-board.
Author: Lori Garzio (lgarzio@marine.rutgers.edu)

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