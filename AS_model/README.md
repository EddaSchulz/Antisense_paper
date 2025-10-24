# Antisense Model

## Description
This folder contains all code that was used to perform and analyze the simulations of the antisense model in Mutzel et al., 2025 

1. *simulations* contains the FST files of the full model and all model simplifications, as well as the bursting and mitosis simulations (under revision).
2. *scripts* contains all code that was used to generate parameter sets, perform and summarize the simulations
3. *figures* contains the antisense model related figures of the paper 

## Information on the simulations
- Scripts to generate the parameter sets are written in matlab
- Scripts to analyze the simulations and generate the figures are written in matlab
- The stochastic antisense model is written in matlab and c++. To allow matlab to call the c++ functions a C++ MEX file of the .cpp file needs to be constructed with:
```
mex reaction_general_AS_model_200319_2.cpp
```
Compilation requires the Boost C++ libraries.


