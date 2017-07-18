# DiagMC
Diagrammatic Monte Carlo

Here is a project to implement the Diagrammatic Monte Carlo method to spin system. First we map spin operators to fermionic ones, then the spin Hamiltonian reduces to dispersionless fermionic Hamiltonian with only interactions. Then the fermion system can be studied by Green's functions, and there comes the Feynman diagrams.

They are several files in this project.

"1storder.cpp": This file munipulate 1st order diagrams. <br />
"calculation.cpp": calculate the relating varibles, such as Green's function, self-energy etc.<br />
"create.cc": create/delete new element in the diagram, such as to create/delete lines.<br />
"fft.cpp": relating to fast Fourier transformations.<br />
"linkedlist.cpp": the functions manipulate the linked-list which construct the diagrams.<br />
"momenta.c": munipulate the momenta on lines.<br />
"newh.h": h file.<br />
"newmain.cpp": main file, run iterations.<br />
"profile.cpp": profile stores data.<br />
"test.c": test the diagram, connection, momentum conservation etc.<br />
"updates.cpp": old version, update the diagrams topology or internal variables.<br />
"updates1.cpp": new version, update the diagrams topology or internal variables.<br />
