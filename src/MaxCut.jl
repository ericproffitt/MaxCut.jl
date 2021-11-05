module MaxCut

using LinearAlgebra
using Convex
using SCS

"The Goemans-Williamson algorithm for the MAXCUT problem."

function maxcut(W::Matrix{<:Real}; iter::Int=100, tol::Real=0)
	"Partition a graph into two disjoint sets such that the sum of the edge weights
	which cross the partition is as large as possible (known to be NP-hard)."

	"A cut of a graph can be produced by assigning either 1 or -1 to each vertex. The Goemans-Williamson 
	algorithm relaxes this binary condition to allow for vector assignments drawn from the (n-1)-sphere
	(choosing an n-1 dimensional space will ensure seperability). This relaxation can then be written as 
	a semidefinite program. Once the optimal vector assignments are found, origin centered hyperplanes are generated and 
	their corresponding cuts evaluated. After 'iter' trials, or when the desired tolerance is reached,
	the hyperplane with the highest corresponding binary cut is used to partition the vertices."
	
	"W:		Adjacency matrix."
	"tol:	Maximum acceptable distance between a cut and the MAXCUT upper bound."
	"iter:	Maximum number of hyperplane iterations before a cut is chosen."

	LinearAlgebra.checksquare(W)
	issymmetric(W)					|| throw(ArgumentError("Adjacency matrix must be symmetric."))
	all(W .>= 0)					|| throw(ArgumentError("Adjacency matrix must be nonnegative."))
	all(iszero.(diag(W)))			|| throw(ArgumentError("Diagonal of adjacency matrix must be zero (no self loops)."))
	(tol >= 0)						|| throw(ArgumentError("The tolerance must be nonnegative."))
	(iter > 0)						|| throw(ArgumentError("The number of iterations must be a positive integer."))

	"This is the standard SDP Relaxation of the MAXCUT problem, a reference can be found at,
	http://www.sfu.ca/~mdevos/notes/semidef/GW.pdf"
	k = size(W, 1)
	S = Semidefinite(k)
	
	expr = dot(W, S)
	constr = [S[i,i] == 1.0 for i in 1:k]
	problem = minimize(expr, constr...)
	solve!(problem, SCS.Optimizer(verbose=false))

	## ensure symmetric positive-definite
	A = 0.5 * (S.value + S.value')
	A += (max(0, -eigmin(A)) + 1e-14) * Matrix(I, size(A, 1), size(A, 1))

	X = Matrix(cholesky(A))

	## a non-trivial upper bound on MAXCUT
	upperbound = (sum(W) - dot(W, S.value)) / 4

	## random origin-centered hyperplanes, generated to produce partitions of the graph
	max_cut = 0
	max_partition = nothing

	for i in 1:iter
		gweval = X' * randn(k)
		partition = (findall(gweval .>= 0), findall(gweval .< 0))
		cut = sum(W[partition...])

		if cut > max_cut
			max_partition = partition
			max_cut = cut
		end

		(upperbound - max_cut < tol) && break
		(i == iter) && println("Max iterations reached.")
	end
	return max_cut, max_partition
end

export maxcut

end
