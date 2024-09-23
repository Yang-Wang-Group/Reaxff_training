# An Improved Reactive Force Field Parameter Optimization Framework Based on Simulated Annealing and Particle Swarm Optimization Algorithms

This package combined the Simulated Annealing algorithm and Particle Swarm Optimization algorithm to optimize the Reactive Force Field, besides, we introduced a concentrated attention mechanism to 
improve accuracy.
This work is completed using C language.

## Usage
In the linux system, after downloading all files, you can use the ' make ' command to compile.

### Description of each file

#### file-XXX-XXX
This file contains the reference density functional theory(DFT) calculation results and the corresponding model information.

#### maxdp.reax
This file contains parameters that needs to be optimized and their maximum step size. If a parameter needs to be optimized, its step size should be set to a positive number, otherwise set to zero.

#### uplimit.reax
This file contains the uplimit of parameters.

#### downlimit.reax
This file contains the uplimit of parameters.

## LAMMPS
LAMMPS needs to be installed in advance

## Contact and questions
If you have any questions, pleas contanct with us at yang.wang@swjtu.edu.cn.
