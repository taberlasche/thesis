# thesis

+ Doc/tex: tex docs for the thesis and markdown docs for drafts / templates
+ Code: implementations
	+ sdp.py: solving the semidefinite program
	+ gram.py: finding the gram vectors
	+ algo.py: randomized rounding algorithm
	+ whole.py: whole program
	+ operations.py: some common operations

Funktionenüberblick:
	+ buildC(n) baut C für einen Hamiltonian der Art X_1X_2+Z_1Z_2+Z_{n+1}X_3+Z_{n+1}X_4...
	+ sdp1 baut die (versuchte) Lösung für das Sdp dieses Hamiltonians
	+ chainC(n) baut C für einen Hamiltonian der Art X_i X_{i+1}, eine Spinkette ohne lineare Terme
	+ sdp baut die (versuchte) Lösung für das SDP dieses Hamiltonians
