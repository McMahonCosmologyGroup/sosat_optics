# Simons Obervatory Small Aperture Telescope (sosat) Optical Simulation

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)

Optical simulation of the Simons Observatory Small Aperture Telescope. <br />
Author: Grace Chesmore

## Installation
- ```git clone git@github.com:McMahonCosmologyGroup/sosat_optics.git```
- ```cd sosat_optics```
- ```pip3 install .```

## Usage
Two notebooks are provided to demonstrate the use of the optical simulation. <br />
### Near Field
This [notebook](https://github.com/McMahonCosmologyGroup/sosat_optics/blob/main/notebooks/sat_nearfield.ipynb) demonstrates the use of the optical simulation to create the near-fields in front of the SAT window. <br />

### Far Field
This [notebook](https://github.com/McMahonCosmologyGroup/sosat_optics/blob/main/notebooks/sat_farfield.ipynb) demonstrates the use of the optical simulation to propagate the near field simulated beam into the far field of the telescope.  With the simulated near-fields $b(x,y)$ above the SAT window, the far-field $B(\theta_x,\theta_y)$ is calculated using the relation:

$$ B(\theta_x,\theta_y) = \int_{aperture} b(x,y)e^{i2\pi(x\theta_x + y\theta_y)} dx dy$$

where $\theta_x$ and $\theta_y$ are the polar coordinates of the far-field, $x$ and $y$ are the coordinates of the near-field, and $b(x,y)$ is the measured near-field.  Here, we are integrating over the area of aperture, which is the Stop of the SAT.

It is important to note that in practice, one can only obtain the far-field beams from a near-field measurement if the near-field data includes phase information.

## Contributions
If you have write access to this repository, please:
* create a new branch
* push your changes to that branch
* merge or rebase to get in sync with main
* submit a pull request on github
* If you do not have write access, create a fork of this repository and proceed as described above. For more details, see Contributing.
