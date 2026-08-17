# StochasticBC
C++ implementation of Stochastic Barycentric Coordinates [de Goes and Desbrun, 2024]


# Update 08/15

Created new repo, inheriting a lot of data structures from Profile Mover, specifically the mesh and curvenet structures.
Wrote some skeleton code in the cagedeformer class; mostly need to update other data structures to accommodate its new queries (bases, closest points, smoothing operators, etc.). Biggest TODO's are to get closest points on curvenets (needs its own BVH and specialized structure for querying closest points + bases), as well as figuring out how to generalize the code such that it works with various different input formats for cage type without needing too many specialized inherited classes.

Surprisingly simple implementation of harmonic coordinates... suspiciously simple, even.