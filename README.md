**The Goemans-Williamson algorithm for the MAXCUT problem.**

Partition a graph into two disjoint sets such that the sum of the edge weights which cross the partition is as large as possible (known to be NP-hard).

A cut of a graph can be produced by assigning either 1 or -1 to each vertex. The Goemans-Williamson algorithm relaxes this binary condition to allow for vector assignments drawn from the (n-1)-sphere (choosing an n-1 dimensional space will ensure seperability). This relaxation can then be written as an SDP. Once the optimal vector assignments are found, origin centered hyperplanes are generated and their corresponding cuts evaluated. After 'iter' trials, or when the desired tolerance is reached, the hyperplane with the highest corresponding binary cut is used to partition the vertices.

### Dependencies
```julia
LinearAlgebra
Convex
SCS
```

### Arguments
```julia
W:		Adjacency matrix.
tol:	Maximum acceptable distance between a cut and the MAXCUT upper bound.
iter:	Maximum number of hyperplane iterations before a cut is chosen.
```

### Example
```julia
	W = [0 5 2 1 0; 
		 5 0 3 2 0; 
		 2 3 0 0 0; 
		 1 2 0 0 4; 
		 0 0 0 4 0]

	max_cut, max_partition = maxcut(W)
	
	@show max_cut
	### max_cut = 14.0

	@show max_partition
	### max_partition = ([1, 3, 4], [2, 5])
```