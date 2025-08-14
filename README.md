## Improving quantum communication rates with permutation-invariant codes

This repository contains MATLAB code used to compute the coherent information of certain permutation-invariant quantum codes using representation theory as described in the paper 

> Sujeet Bhalerao, Felix Leditzky: "Improving quantum communication rates with permutation-invariant codes", [arXiv:2508.09978](https://arxiv.org/abs/2508.09978).

### Requirements

- Parallel Computing Toolbox
- MATLAB's particleswarm function, available in the "Global Optimization Toolbox"
- `fastexpm` function for computing matrix exponentials of sparse matrices (available at MATLAB File Exchange [here](https://www.mathworks.com/matlabcentral/fileexchange/84058-fast-exponential-matrix-for-matlab-full-sparse-fastexpm))



### Usage

#### To optimize over permutation-invariant codes:

1. Run the optimization script `pauli_opt.m` with the desired parameters. This is currently set up to optimize over permutation-invariant codes for the $2$-Pauli channel with $9$ copies at a noise level of $0.2271$, which is slightly beyond the hashing bound for this channel.

2. The results will be saved in the `optimization_results` directory.


#### To compute coherent information for a specific permutation-invariant code:

Run the script `single_code_ci.m` with the desired parameters and input density matrices. The current parameters are set for $9$ copies of a $2$-Pauli channel using the code listed in Table 4 of the paper. 

- Using purifications:
    For qubit channels, this can compute the coherent information for a single code for up to $100$ channel qubits. For channels with qutrit output (e.g. the dephrasure channel), this can compute the coherent information for a single code for up to $30$ channel qubits.

- Using complementary channels: 
    For channels with a $3$-dimensional environment (e.g. the $2$-Pauli channel), this can compute the coherent information for a single code for upto $40$ channel qubits. For channels with a $4$-dimensional environment (e.g. the BB84 channel), this can compute the coherent information for a single code for upto $15$ channel qubits.


### Files

#### Optimization

- `pauli_opt.m`: Main optimization script that sets up the problem and calls the particle swarm optimization function.

- `compute_ci_symmetries_opt.m`: Objective function for the particle swarm optimization, which computes the coherent information for a given set of input states and probabilities using the complementary channel.

- `compute_entropy_d2.m`, `compute_entropy_d3.m`, `compute_entropy_d4.m`: Functions to compute the entropy of states in dimensions $2$, $3$, and $4$ respectively (e.g. as outputs of qubit channels with environments of dimensions  $3$ and $4$) using representation theory.

#### Evaluating specific codes
- `eval_codes.m`: Script to verify the coherent information for the specific codes listed in Table 4 of the paper. Each case in the script sets the required channel parameters and code states, then calls the appropriate computation function. This script can be used as a template to evaluate other codes by modifying the parameters within a 'case' block.

- `compute_ci_symmetries.m`: Function to compute coherent information for a given set of input states and probabilities using our representation-theoretic approach - makes use of the complementary channel.
- `compute_ci_purification.m`: Function to compute coherent information using purifications for a given set of pure input states and probabilities using our representation-theoretic approach. 
- `compute_ci_purification_qutrit.m`: Similar to the function above but for channels with qutrit output.

#### Auxiliary Files

- `compute_symm_rep_formula.m`: computes matrix elements for the GL(2) irrep on the symmetric subspace (see Appendix B in the paper).

- `compute_S_RA_GL2.m`, `compute_S_RA_log.m`, `compute_S_RA_purification_qutrit.m`: Functions to compute entropy at the environment for the coherent information. 
`compute_S_RA_GL2` uses explicit GL(2) irreps, `compute_S_RA_log` uses a log/exp approach, and `compute_S_RA_purification_qutrit` is for qutrit channels.

- `get_states_bloch.m` and `get_states_measure.m`: two parametrization schemes for qubit input states, one using Bloch vectors and the other using measurement outcomes. See Appendix D of the paper for details.

- `comp_kraus.m`: Function to compute the Kraus operators of the complementary channel for a given set of operators and output dimensions.

- `q_logm.m`: Function to compute logarithm of the positive part of a matrix.

- The folders `E_matrices` and `E_matrices_sparse` contain precomputed matrices for the representation of the Lie algebra $\mathfrak{gl}_d(\mathbb{C})$ as described in the paper. 

- The folders `dimWlambda` and `gt_patterns` contain scripts for computing the dimensions of Specht modules and generating Gelfand-Tsetlin patterns for a given partition, respectively.

#### SageMath Scripts

To install SageMath, follow the instructions at [SageMath Installation](https://doc.sagemath.org/html/en/installation/index.html).

- `dimWlambda/spechtmodule_dim.py`: Script to compute the dimension of the Specht module for a given partition.
- `gt_patterns/gen_gt_patterns.py`: Script to generate Gelfand-Tsetlin patterns for a given partition and number of rows. 

