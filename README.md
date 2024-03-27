![MuSSIC's Logo](logo/MuSSIC_logo.png)

# MuSSIC
-----------------------
Description of the Code
-----------------------
MuSSIC is a computational tool  for calculation of neutron scattering for coarse grain (CG) systems.
It is based on the method described by Alan Soper in his journal article on coarse-grained Empirical Potential Structure Refinement
(Biochimica et Biophysica Acta 1861, 6,2017,Pages 1652-1660).
The code has been tested for different coarse grain models on surfactants in water to complex polymer melt systems(examples included in usability_tests/). 
Further usability tests on different systems are being conducted to test the generality and applicability of the method to other systems.
The current version of the code (v1.0) can be  used on polymer metls or surfactants forming micelles (like SDS, PEG and CTAB). 
Any development  or improvement in the code/method will be updated on Github (link to the page).
The code computes the partial pair radial distribution functions and the partial pair structure factors for the trajectory provided in xyz file format. 
The total struture factor F_CG(Q) is the weighted sum of all  partial pair structure factors(detailed notes is given in the documentation).
The userguide 'MuSSIC_user_guide' provide details on compilation, input files and output files.

-----
# User guide
-----
Detailed user guide together with example input files is available [here]([https://github.com/disorderedmaterials/MuSSIC/blob/master/MuSSIC_user_guide_v1.1.pdf])

MuSSIC is wriiten in FORTRAN and requires gfortran compiler with openmp support for parallel implementation. 

Compile the code by typing 'make' in the sourcecode and copy the executable 'scattcg' to the desired folder to run the simulations

>make

>cp scattcg /desired_location

>cd /desired_location

>./scattcg

-------------
# Documentation
-------------

Documentation on the neutron scattering calculation for coarse grained simulations and validation tests is available [here](https://github.com/disorderedmaterials/MuSSIC/blob/master/MuSSIC_documentation.pdf).

Documentation on the tests from coarse-grained simulations is a work in progess. Details will be added as the work progresses.

