![MuSSIC's Logo](icon/MuSSIC_logo.png)

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
usage
-----
Detailed user guide together with comprehensive examples is available at 
To compile the code, type 'make'
To run the code, type './scattcg'

-----------
INPUT files
-----------
input.dat    -  main input file (details explained below).
model.dat    -  model parameters of the system.
CG_trajectory   -  coarse-grained trajectory in xyz format 

------------
OUTPUT files
------------
#File name          #description

fort.7               # log file
unweighted-*-*.gr    # partial pair g_{ij}(r) (* - atom types in the pair)
unweighted-*-*.sq    # partial pair S_{ij}(Q) (* -atom types in the pair)
unweigthed-total.gr  #  unweighted S(Q)
unweighted-total.sq  # unweighted G(r)
**-weighted-*-*.sq   # weighted partial pair S_{ij}(Q)(** -isotpologue name, * - atom types in the pair)
**-weighted-*-*.gr   # weighted partial pair g_{ij}(r)(** -isotpologue name, * - atom types in the pair)
**-weighted-total.sq # weighted S(Q) (**-isotopologue name)
**-weighted-total.gr # weighted G(r) (**-isotopologue name)

----------------------------
Main input  file 'input.dat'
----------------------------
 
title:            # title of the system (optional)
Methanol-water
ntypes            # number of atom types
6
nframes:          # number of frames in the trajectory
5
box_length:       # sides of the box in Angstroms
35.5689d0  35.5689d0  35.5689d0
bin_width:        # bin width for computing g(r) in Angstroms
0.03d0
cutoff:           # half of the box length snapped to bin width
17.759d0
label_atom:       # atom label given in the trajectory
OW
HW
CT
HC
OH
HO
scattering_lengths: # bound coherent scattering lenghts for the calculation of neutron scattering
5.803d0
-3.7390d0
6.646d0
-3.739d0
5.803d0
-3.739d0
6.671d0
exchangeable_atoms:
2                    # number of exchangeable atoms
HW                   # atoms labels of exchangeables
HO
number_isotopologues:
1
isotopologues:       # description of isotope ratio for all atom types in the isotopologue (not relevant but provide  1.0 for all  natural and 0.0 to isotope)
sample_1
atom    natural   isotope
OW       1.0      0.0
HW       0.0      1.0
CT       1.0      0.0
HC       1.0      0.0
OH       1.0      0.0
HO       1.0      0.0
Isotopic_substitution:
manual                # isotopic substitution can be implemented either by choosing                    
                              ‘manual’ or ‘ratios’
                           #  ‘manual’ option is be used to generate manually deuterated trajectories. Test case 4 shown in the documentation explains the results with manual deuteration for PA66. While the test case 5 shows the results obtained using ‘ratios’ option for C10TAB in water (Refer to verification tests in the documentation)
 (note: the current version of the code supports only one isotope of the atoms and only one isotopologue can be tested per run if using manual substitution)

qmin:                # Q minimum for the calculation S(Q)
0.1d0
qdelta:              # bin width for Q in (1/Angstroms)
0.05d0
qmax:                # Q maximum
30.0d0
broadening:             # broadening applied for intra molecular correlations
OmegaDependentGaussian
FWHM:
0.02d0
windowfn:               # Lorch function to reduce the termination ripples while performing back Fourier transform S(Q) -> G(r)
Lorch0
endinput

----------------------
Model file 'model.dat'
----------------------
num_species   # number of molecule types
2
nmol_species  # number of molecules of each type
300
600
size_per_species # numbe rof atoms in each molecule
3
6
sequence_atoms  # sequence of atoms in molecules eg: 3 atoms in 1st molecule and 6 atoms in second molecule
OW
HW
HW
CT
HC
OH
HC
HC
HO