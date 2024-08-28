# Tuning higher order Van Hove singularities in Strontium Ruthenate (Sr<sub>2</sub>RuO<sub>4</sub>)
The folder **Sr2RuO4** contains the Mathematica notebooks and related material for https://arxiv.org/abs/2310.15331
Most of the code is written in Wolfram Mathematica 13, with some parts written C++ (requires LAPACK++ and OpenMP).
The Wannierization of the DFT calculations was done by Luke Rhodes. The code requires my Mathematica Package, **BandUtilities**. This is contained in the folder **Package** along with some illustrative examples located in the folder **Examples**. 

The notebooks in **Sr2RuO4** are organised by the type of calculations they contain, with the title reflecting that:
- 1_ARPES_best_fit.nb
In this notebook, the Wannierised tight-binding models are loaded and polynomial interpolated with symmetrisation with the help of the functions defined in BandUtilities. This generates a ϴ dependant tuneable model that is then fit to ARPES data to obtain 1) overall band renormalisation, 2) Spin orbit strength and 3) Fermi level adjustment.
- 
Please reach out if you have any questions.

Anirudh Chandrasekaran

5 August 2024
