## Global Optimization for Cardinality-constrained Minimum Sum-of-Squares Clustering via Semidefinite Programming

<p align="center">
  <img src="https://github.com/antoniosudoso/cc-sos-sdp/blob/main/logo.svg" width="200" height="200" />
</p>


**CC-SOS-SDP** is an exact algorithm based on the branch-and-cut technique for solving the Minimum Sum-of-Squares Clustering (MSSC) problem with cardinality constraints described in the paper ["Global Optimization for Cardinality-constrained Minimum Sum-of-Square Clustering via Semidefinite Programming"](https://doi.org/10.1007/s10107-023-02021-8). This repository contains the C++ source code, the MATLAB scripts, and the datasets used for the experiments.

> Piccialli, V., Sudoso, A.M. Global optimization for cardinality-constrained minimum sum-of-squares clustering via semidefinite programming. Math. Program. (2023). https://doi.org/10.1007/s10107-023-02021-8

## Installation
**CC-SOS-SDP** calls the semidefinite programming solver [SDPNAL+](https://blog.nus.edu.sg/mattohkc/softwares/sdpnalplus/) by using the [MATLAB Engine API](https://www.mathworks.com/help/matlab/calling-matlab-engine-from-cpp-programs.html) for C++. It requires the MATLAB engine library *libMatlabEngine* and the Matlab Data Array library *libMatlabDataArray*. **CC-SOS-SDP** calls the integer programming solver [Gurobi](https://www.gurobi.com/). **CC-SOS-SDP** uses [Armadillo](http://arma.sourceforge.net/) to handle matrices and linear algebra operations efficiently. Before installing Armadillo, first install OpenBLAS and LAPACK along with the corresponding development files. **CC-SOS-SDP** implements a configurable thread pool of POSIX threads to speed up the branch-and-bound search.

Ubuntu and Debian instructions:
1) Install MATLAB (>= 2016b)
2) Install Gurobi (>= 9.1)
3) Install CMake, OpenBLAS, LAPACK and Armadillo:
 ```
sudo apt-get update
sudo apt-get install cmake libopenblas-dev liblapack-dev libarmadillo-dev
```
4) Open the makefile `clustering_c++/Makefile` 
	- Set the variable `matlab_path` with your MATLAB folder.
5) Compile the code:

```
cd clustering_c++/
make
```

4) Download SDPNAL+, move the folder `clustering_matlab` containing the MATLAB source code of **CC-SOS-SDP** in the SDPNAL+ main directory and set the parameter `SDP_SOLVER_FOLDER` of the configuration file accordingly. This folder and its subfolders will be automatically added to the MATLAB search path when **CC-SOS-SDP** starts.

The code has been tested on Ubuntu Server 20.04 with MATLAB R2020b, Gurobi 9.5 and Armadillo 10.2.

## Configuration
Various parameters used in **CC-SOS-SDP** can be modified in the configuration file `clustering_c++/config.txt`:

- `BRANCH_AND_BOUND_TOL` - optimality tolerance of the branch-and-bound
- `BRANCH_AND_BOUND_PARALLEL` -  thread pool size: single thread (1), multi-thread (> 1)
- `BRANCH_AND_BOUND_MAX_NODES` - maximum number of nodes
- `BRANCH_AND_BOUND_VISITING_STRATEGY` - best first (0),  depth first (1), breadth first (2)
- `MATLAB_SESSION_THREADS_ROOT` - number of threads for the MATLAB session at the root
- `MATLAB_SESSION_THREADS_CHILD` - number of threads for the MATLAB session of children nodes
- `SDP_SOLVER_FOLDER` - full path of SDPNAL+ folder
- `SDP_RELAXATION` - vector lifting SDP (0), matrix lifting SDP (1)
- `SDP_SOLVER_TOL` - accuracy of SDPNAL+
- `SDP_SOLVER_MAX_ITER` - maximum number of iterations
- `SDP_SOLVER_MAX_TIME` - maximum time in seconds
- `SDP_SOLVER_VERBOSE` - do not display log (0), display log (1)
- `CP_MAX_ITER_ROOT` - maximum number of cutting-plane iteration root node
- `CP_MAX_ITER_CHILD` -  maximum number of cutting-plane iteration child node
- `CP_TOL` - tolerance between two consecutive cutting-plane iterations
- `CP_MAX_INEQ` - maximum number of valid inequalities to separate
- `CP_PERC_INEQ` - fraction of the most violated inequalities to add
- `CP_INHERIT_INEQ` - do not inherit inequalities (0), inherit inequalities (1)
- `CP_EPS_INEQ` - tolerance for checking the violation of the inequalities
- `CP_EPS_ACTIVE` - tolerance for detecting active inequalities
- `GUROBI_FOLDER` - gurobi solver path
- `HEURISTIC_VERBOSE` - do not display log (0), display log (1)

## Usage
```
cd clustering_c++/
./bb <DATASET> <K> <C_1> <C_2> ... <C_K> <LOG> <RESULT>
```
- `DATASET` - path of the dataset
- `K` - number of clusters
- `C_1 C_2 ... C_K` - cluster sizes (cardinality constraints)
- `LOG` - path of the log file
- `RESULT` - path of the optimal cluster assignment matrix

File `DATASET` contains the data points `x_ij` and the must include an header line with the problem size `n` and the dimension `d`:

```
n d
x_11 x_12 ... x_1d
x_21 x_22 ... x_2d
...
...
x_n1 x_n2 ... x_nd
```

## Log

The log file reports the progress of the algorithm:

- `N` - size of the current node
- `NODE_PAR` - id of the parent node
- `NODE` - id of the current node
- `LB_PAR` - lower bound of the parent node
- `LB` - lower bound of the current node
- `FLAG` - termination flag of SDPNAL+
    -  `0` - SDP is solved to the required accuracy
    -  `1` - SDP is not solved successfully
    -  `-1, -2, -3` - SDP is partially solved successfully
- `TIME (s)` - running time in seconds of the current node
- `CP_ITER` - number of cutting-plane iterations
- `CP_FLAG` - termination flag of the cutting-plane procedure
    - `-3` - current bound is worse than the previous one
    - `-2` - SDP is not solved successfully
    - `-1` - maximum number of iterations
    -  `0` - no violated inequalities
    -  `1` - maximum number of inequalities
    -  `2` - node must be pruned
    -  `3` - cutting-plane tolerance
- `CP_INEQ` - number of inequalities added in the last cutting-plane iteration
- `PAIR TRIANGLE CLIQUE` - average number of added cuts for each class of inequalities
- `UB` - current upper bound
- `GUB` - global upper bound
- `I J` - current branching decision
- `NODE_GAP` - gap at the current node
- `GAP` - overall gap 
- `OPEN` - number of open nodes

Log file example:

```
DATA_PATH: /home/ubuntu/CC-SOS-SDP/instances/iris.txt 150 4 3
CLUSTER_SIZES: 50 50 50 
LOG_PATH: /home/ubuntu/CC-SOS-SDP/log_vl/log_vl_iris.txt

BRANCH_AND_BOUND_TOL: 0.0001
BRANCH_AND_BOUND_PARALLEL: 2
BRANCH_AND_BOUND_MAX_NODES: 200
BRANCH_AND_BOUND_VISITING_STRATEGY: 0

MATLAB_SESSION_THREADS_ROOT: 14
MATLAB_SESSION_THREADS_CHILD: 7

SDP_SOLVER_FOLDER: /home/ubuntu/CC-SOS-SDP/SDPNAL+/
SDP_RELAXATION: 0 (VECTOR LIFTING)
SDP_SOLVER_TOL: 0.0001
SDP_SOLVER_MAX_ITER: 50000
SDP_SOLVER_MAX_TIME: 7200
SDP_SOLVER_VERBOSE: 0

CP_MAX_ITER_ROOT: 20
CP_MAX_ITER_CHILD: 10
CP_TOL: 0.0001
CP_MAX_INEQ: 100000
CP_PERC_INEQ: 0.1
CP_INHERIT_INEQ: 1
CP_EPS_INEQ: 0.0001
CP_EPS_ACTIVE: 0.0001

GUROBI_FOLDER: /home/ubuntu/gurobi952/
HEURISTIC_VERBOSE: 0


|    N| NODE_PAR|    NODE|      LB_PAR|          LB|    FLAG|  TIME (s)| CP_ITER| CP_FLAG|   CP_INEQ|          UB|         GUB|     I      J|     NODE_GAP|          GAP|  OPEN|
|  150|       -1|       0|        -inf|     81.2778|      -1|         4|       0|       1|         0|     81.2778|    81.2778*|    -1     -1|  1.05299e-09|  1.05299e-09|     0|

TIME: 5 sec
NODES: 1
ROOT_GAP: 1.05299e-09
GAP: 0
OPT: 81.2778

```
