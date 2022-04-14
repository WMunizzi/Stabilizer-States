# StabilizerStateData
Data files accompanying paper (ArXiv #) by Cynthia Keeler, Jason Pollack, and William Munizzi.

StabilizerState file contains the full set of stabilizer states for 1,2,3,4, and 5 qubits. Each set is contained in its own .wl file, able to be imported directly into Mathematica notebook. States containined in each file are organized as length 2^n bit-addresses. A bit-address is interpreted as the ordered set of coefficients multiplying each basis ket of an n-qubit system, e.g. the bit-address (1,0,0,1,0,0,i,i) indicates the state \ket{000}+\ket{011}+i\ket{110}+i\ket{111}, a 3-qubit state. Note: The ordering of qubits within a ket is given from RIGHT TO LEFT, e.g. the rightmost digit corresponds to the first qubit of the system, while the leftmost digit represents the n^th qubit of an n-qubit syste. The binary string within each ket directly equates the left-acting tensor product of single-qubit states (0 and 1).



File repository consisting of stabilizer states generated to accompany (paper). The full set of stabilizer states is included up through 5 qubits.
