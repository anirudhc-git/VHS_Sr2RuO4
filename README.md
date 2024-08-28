# Tuning higher order Van Hove singularities in Strontium Ruthenate (Sr<sub>2</sub>RuO<sub>4</sub>)
The folder **Sr2RuO4** contains the Mathematica notebooks and related material for https://arxiv.org/abs/2310.15331
Most of the code is written in Wolfram Mathematica 13, with some parts written C++ (requires LAPACK++ and OpenMP).
The Wannierization of the DFT calculations was done by Luke Rhodes. I use a wrapper for LAPACK's zheev that was written by R Ganesh, based on an earlier version by Scott Shaw. The Mathematica notebooks require my Mathematica Package, **BandUtilities**. This is contained in the folder **Package** along with some illustrative examples located in the folder **Examples**. 

The notebooks in **Sr2RuO4** are organised by the type of calculations they contain, with the title reflecting that:
- 1_ARPES_best_fit.nb: 
In this notebook, the Wannierised tight-binding models are loaded and polynomial interpolated with symmetrisation with the help of the functions defined in BandUtilities. This generates a ϴ dependent tuneable model that is then fit to ARPES data to obtain 1) overall band renormalisation, 2) Spin orbit strength and 3) Fermi level adjustment.
- 2_HOVHS.nb:
Using the ϴ dependent model, we examine the tuning and evolution of the higher order Van Hove singularity with ϴ and nematicity and generate the phase diagrams used in the main text.
- 3_M_point_staggered_chemical_potential.nb: 
The addition of staggered chemical potential instead of nematicity, leads to the formation of a four fold Higher order Van Hove singularity. This situation is analysed here.
- 4_Log_to_Power_law.nb:
The transition from logarithmic to power law density of states in the vicinity of the critical tuning (to a higher order Van Hove singularity) is worked out in this notebook.

The notebooks are to be executed in order. They are amply commented and contain step by step instructions to aid the user. Please reach out if you have any questions.

Anirudh Chandrasekaran

5 August 2024
