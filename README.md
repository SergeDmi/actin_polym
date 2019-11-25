# Actin polymerization simulation
A simple code to simulate actin monomer diffusion and reaction with filaments.

## Installation

You need a working python environment. Then install pandas, numpy, and yaml :  
$ pip3 install numpy --user
$ pip3 install yaml --user
$ pip3 install pandas --user

Then download actin_polymerization.py from this github repo.

## Usage

python3 actin_polymerization.py config.yaml

## Configuration files

config.yaml is a suitable yaml-type config file, e.g.  
```yaml
simulation:
    Tmax: 1000
    N_monomers: 100000
    box: [150,120,120]
    success_frac: 1.0
filaments:
    f0:
        position: [60,19]
        orientation: [1,0]
        length: 50
    f1:
        position: [60,20]
        orientation: [1,0]
        length: 100
```

# Serge Dmitrieff -- http://biophysics.fr
