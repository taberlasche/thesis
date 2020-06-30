# TODO
+ [ ] states
+ [x] positive semidefinite matrices properties
+ [ ] algorithms
+ [ ] pauli matrices
+ [ ] gellmann matrices
+ [ ] qudits
+ [ ] bloch sphere qubits
+ [ ] bloch sphrer qutrits
+ [ ] non positiv-semidefinite states

# States
+ pos. semidefinite
+ trace 1
+ hermitian

# Properties positive semidefite matrices $M$

$M$ has non-zero eigenvalues

The diagonal entries of of a positive definite matrix are real and non-negative
An $n\times n$ hermitian matrix  is pos. def. if
  $tr(M) > 0$
 and $\frac{(tr(M))^2}{tr(tr(M^2)} > n-1$

$M,N\ge 0 \implies M\otimes N \ge 0$

A hermitian matrix $M$ is positive definite if and only if it has a unique Cholesky decomposition, i.e. the matrix $M$ is positive definite if and only if there exists a unique lower triangular matrix $L$, with real and strictly positive diagonal elements, such that M = LL^{*)
Furthermore, positive-semidefinite Matrices have such a decomposition if the diagonal elements of $L$ are allowed to be zero.



# Pauli Matrices
 $$
 X = \begin{bmatrix} 0 & 1 \\
                    1 & 0
        \end{bmatrix},~

    ~Y = \begin{bmatrix} 0 & -i \\
                    i & 0
         \end{bmatrix},~

    ~Z = \begin{bmatrix} 1 & 0 \\
                    0 & -1
        \end{bmatrix}
$$
Properties:
+ traceless
+ unitary
+ hermitian
