# boarbitrary
Efficient Framework for Solving Plasma Waves Kinetic Dispersion Relation with Arbitrary Distributions

This ('BO-Arbitrary') is an extension of the kinetic electromagnetic magnetized dispersion relation solver PDRK/BO (https://github.com/hsxie/pdrk, https://github.com/hsxie/bo) to arbitrary distributions, with three versions, Hermite-Hermite (HH, J-pole matrix similar to Maxwellian case), GPDF-Hermite (GH, more accurate for parallel integral), GPDF-GPDF (GG, can use FFT).

fsolve root finding versions are provided for all the above three cases, i.e., boem3dHHroot, boem3dGHroot and boem3dGGroot. 

Matrix version to obtain all solutions is provided for HH case.

Ref: H. S. Xie, Efficient Framework for Solving Plasma Waves with Arbitrary Distributions, 2025, https://arxiv.org/abs/2501.06477. Also, Plasmas 32, 060702 (2025); [doi: 10.1063/5.0275307](https://doi.org/10.1063/5.0275307).

8:41 2025/1/10


Julia version for BO-arbitrary Hermite-Hermite (HH) version, and product bi-kappa (BO-PBK) version can be found at: https://juliaspacephysics.github.io/PlasmaBO.jl/dev (By Zijin ZHANG at UCLA). 

Matlab product bi-kappa (BO-PBK) version can be found at: https://github.com/baiweiphys/BOPBK (By Wei BAI).

7:23 2025/12/27
