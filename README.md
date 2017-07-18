# DiagMC
Diagrammatic Monte Carlo

Here is a project to implement the Diagrammatic Monte Carlo method to spin system. First we map spin operators to fermionic ones, then the spin Hamiltonian reduces to dispersionless fermionic Hamiltonian with only interactions. Then the fermion system can be studied by Green's functions, and there comes the Feynman diagrams.

They are several files in this project.

"1storder.cpp": This file munipulate 1st order diagrams._
"calculation.cpp": calculate the relating varibles, such as Green's function, self-energy etc. 

"create.cc": create/delete new element in the diagram, such as to create/delete lines.

"fft.cpp": relating to fast Fourier transformations.

"linkedlist.cpp": the functions manipulate the linked-list which construct the diagrams.

"momenta.c": munipulate the momenta on lines.

"newh.h": h file.
"newmain.cpp": main file, run iterations.
"profile.cpp": profile stores data.
"test.c": test the diagram, connection, momentum conservation etc.
"updates.cpp": old version, update the diagrams topology or internal variables.
"updates1.cpp": new version, update the diagrams topology or internal variables.
