# pppmodel

pppmodel is a reimplementation of a model orginally proposed by McIntyre et al. (1989) describing a F-type non-oxidative pentose phosphate pathway in erythrocytes. Later this model was adapted by Berthon et al. (1993) for the in silico replication of NMR studies.

Install with pip install -e .

This repository/ package is an application example of the Python package "modelbase" for simulating ordinary differential equations in metabolic networks. It is shown that the "modelbase" is capable of inferring the dynamics of labels in metabolic pathways (here the non-oxidative pentose phosphate pathway).

The folder pppmodel contains the basic model of the non-oxidative pentose phosphate pathway (PPPmodel.py) originally proposed by McIntyre et al. (1989). This model can be simulated using the methods described in the file simulate.py.

However, the reimplementation of this basic model and its expansion to a system that is able to simulate the dynamics of labels using the Python package "modelbase" is shown in the file "LabelPPPmodel.py".

In the folder "examples" you can find a file, "Berthon1993.py", for the reproduction of some figures of Berthon et al. (1993).

Please execute the file "Berthon1993.py" after  installing the Python package modelbase and all its dependencies as well as the package contained in this repository (Installation, see above).
 

## References

H. A. Berthon, W. A. Bubb, and P. W. Kuchel. '$^13$'C n.m.r. isotopomer and
computer-simulation studies of the non-oxidative pentose phosphate pathway of
human erythrocytes. Biochemical Journal, 296(2):379–387, Dec. 1993.

L. M. Mcintyre, D. R. Thorburn, W. A. Bubb, and P. W. Kuchel. Comparison of
computer simulations of the F-type and L-type non-oxidative hexose monophos-
phate shunts with 31p-NMR experimental data from human erythrocytes. Euro-
pean Journal of Biochemistry, 180(2):399–420, Mar. 1989.

