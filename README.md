# thesis

+ Doc/tex: tex docs for the thesis and markdown docs for drafts / templates
+ Code: implementations
	+ sdp.py: solving the semidefinite program
	+ gram.py: finding the gram vectors
	+ algo.py: randomized rounding algorithm
	+ whole.py: whole program
	+ operations.py: some common operations

+Funktionenüberblick:
	+ buildC(n) baut C für einen Hamiltonian der Art X_1X_2+Z_1Z_2+Z_{n+1}X_3+Z_{n+1}X_4...
	+ sdp1 baut die (versuchte) Lösung für das Sdp dieses Hamiltonians. (Beispiele in sdp1_examples))
	+ chainC(n) baut C für einen Hamiltonian der Art X_i X_{i+1}, eine Spinkette ohne lineare Terme
	+ sdp2 baut die (versuchte) Lösung für das SDP dieses Hamiltonians. (Beispiele in sdp2_examples)
	+ test(n) löst das jeweilige sdp mit mc(C) und rundet werte <10e-8 auf 0. Beispiele für die jeweiligen Cs in "buildC_sdp_rounded" (analog zu sdp1(n)) und "chainC_sdp_rounded" (analog zu sdp2(n))
