using Convex
using SCS

"MAXCUT";

"Partition a graph into two disjoint sets such that the sum of the edge weights
from all edges which cross the partition is as large as possible (known to be NP-hard).";

function goemansWilliamson{T <: Real}(W::Array{T, 2}; tol::Real=1e-1, iter::Int=100)
	"A cut of a graph can be produced by assigning either 1 or -1 to each vertex.  The Goemans-Williamson 
	algorithm relaxes this binary condition to allow for vector assignments drawn from the (n-1)-sphere
	(choosing an n-1 dimensional space will insure seperability).  This relaxation can then be written as 
	an SDP.  Once the optimal vector assignments are found, origin centered hyperplanes are generated and 
	their corresponding cuts evaluated.  After 'iter' trials, or when the desired tolerance is reached,
	which ever comes first, the hyperplane with the highest corresponding binary cut is used to partition 
	the vertices.";
	"W:	Adjacency matrix.";
	"tol:	Maximum acceptable distance between a cut and the MAXCUT upper bound.";
	"iter:	Maximum number of hyperplane iterations before a cut is chosen.";
	LinAlg.chksquare(W)
	LinAlg.issym(W)		|| error("Adjacency matrix must be symmetric.")
	all(W .>= 0)		|| error("Entries of the adjacency matrix must be nonnegative.")
	all(diag(W) .== 0)	|| error("Diagonal entries of adjacency matrix must be zero.")
	tol > 0			|| error("The tolerance 'tol' must be positive.")
	iter > 0		|| error("The number of iterations 'iter' must be at least one.")

	"This is the standard SDP Relaxation of the MAXCUT problem, a reference can be found at
	http://www.sfu.ca/~mdevos/notes/semidef/GW.pdf.";
	const k = size(W, 1)
	S = Semidefinite(k)
	
	expr = vecdot(W, S)
	constr = [S[i,i] == 1.0 for i in 1:k]
	problem = minimize(expr, constr...)
	solve!(problem, SCSSolver(verbose=0))
	
	X = full(cholfact(S.value, :U, Val{true}))
	upperbound = (sum(W) - vecdot(W, S.value))/4 # A non-trivial upper bound on MAXCUT.

	"Random origin-centered hyperplanes, generated to produce partitions of the graph.";
	maxcut = 0
	maxpartition = nothing

	for i in 1:iter
		eval = X' * randn(k)
		partition = (find(eval .>= 0), find(eval .< 0))
		cut = sum(W[partition...])

		if cut > maxcut
			maxpartition = partition
			maxcut = cut
		end
		
		if (upperbound - maxcut) < tol; break; end
		if i == iter; println("Max iterations reached."); end
	end
	return round(maxcut,3), maxpartition
end

function test()
	W = [0 5 2 1 0; 
	     5 0 3 2 0; 
	     2 3 0 0 0; 
	     1 2 0 0 4; 
	     0 0 0 4 0]
	maxcut, maxpartition = goemansWilliamson(W)
	println(maxcut)
	println(maxpartition)
end



