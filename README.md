# boarbitrary
Efficient Framework for Solving Plasma Waves Kinetic Dispersion Relation with Arbitrary Distributions

This ('BO-Arbitrary') is an extension of the kinetic electromagnetic magnetized dispersion relation solver PDRK/BO (https://github.com/hsxie/pdrk, https://github.com/hsxie/bo) to arbitrary distributions, with three versions, Hermite-Hermite (HH, J-pole matrix similar to Maxwellian case), GPDF-Hermite (GH, more accurate for parallel integral), GPDF-GPDF (GG, can use FFT).

fsolve root finding versions are provided for all the above three cases, i.e., boem3dHHroot, boem3dGHroot and boem3dGGroot. 

Matrix version to obtain all solutions is provided for HH case.

Ref: H. S. Xie, Efficient Framework for Solving Plasma Waves with Arbitrary Distributions, 2025, https://arxiv.org/abs/2501.06477.

8:41 2025/1/10
