**The Goemans-Williamson algorithm for the MAXCUT problem.**

Partition a graph into two disjoint sets such that the sum of the edge weights which cross the partition is as large as possible (known to be NP-hard).

A cut of a graph can be produced by assigning either 1 or -1 to each vertex. The Goemans-Williamson algorithm relaxes this binary condition to allow for vector assignments drawn from the (n-1)-sphere (choosing an n-1 dimensional space will ensure seperability). This relaxation can then be written as a semidefinite program. Once the optimal vector assignments are found, origin centered hyperplanes are generated and their corresponding cuts evaluated. After `iter` trials, or when the desired tolerance is reached, the hyperplane with the highest corresponding binary cut is used to partition the vertices.

### Installation
```julia
(@v1.6) pkg> add https://github.com/ericproffitt/MaxCut.jl
```

### Dependencies
```julia
LinearAlgebra
Convex
SCS
```

### Arguments
```julia
W:      (positional arg) Adjacency matrix.
tol:    (keyword arg) Maximum acceptable distance between a cut and the MAXCUT upper bound (default=0).
iter:   (keyword arg) Maximum number of hyperplane iterations before a cut is chosen (default=100).
```

### Example
```julia
W = [0 5 2 1 0; 
     5 0 3 2 0; 
     2 3 0 0 0; 
     1 2 0 0 4; 
     0 0 0 4 0];

max_cut, max_partition = maxcut(W);
	
@show max_cut;
## max_cut = 14

@show max_partition;
## max_partition = ([1, 3, 4], [2, 5])
```

### Reference
http://www.sfu.ca/~mdevos/notes/semidef/GW.pdf