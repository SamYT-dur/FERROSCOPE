# FERROSCOPE
Repository for the Evans group ferroelectric search program using the CSD Python API. Only works on Linux.

Every polar entry from the CSD is considered as a potential ferroelectric material.\
Structural modifications are made to each. \
A higher-symmetry structure is searched for using FINDSYM.\
Successful entries are stored in an sqlite database.

Article currently in review.

/PSTG_ONLY contains a stripped down version of findsym for those only interested in erroneous space group choice in the CSD

## Installation
* Download and install the CSD Python API. FERROSCOPE was tested on the Sep24 version.
* Use pip to install extra packages: pymatgen, pyxtal, qc-iodata, stopit, rdkit, miniball, gemmi. Installing these should install all other dependecies e.g. numpy/scipy.
* Install ISOTROPY Linux suite to get FINDSYM locally: https://stokes.byu.edu/iso/isolinux.php
* Ensure both CSD python API and ISOTROPY are referenced in bashrc. e.g. "export ISODATA=/home/user/iso/" and "alias ccdcpython='source /home/user/CCDC/ccdc-software/csd-python-api/miniconda/bin/activate'"

## Running
* MASTER.py contains options for running the search including what structural modifications should be performed.
* Command line. To run type "ccdcpython" to activate CSD Python API.
* "python MASTER.py" should begin the search. Options may pop up to the user depending on settings.
 
## code structure
MASTER.py \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--> SEARCH.py \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--> STRUCTURE.py\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--> TOOLS.py \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--> SYMMETRY_DETECTION.py\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--> SUPERIMPOSE.py



