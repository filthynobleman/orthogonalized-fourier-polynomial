This code implements a basic version of the algorithm described in the paper:

"Informative Descriptor Preservation via Commutativity for Shape Matching,"
Dorian Nogneng and Maks Ovsjanikov, Proc. Eurographics 2017

To try it, simply run test_compute_faust.m in MATLAB. This should produce
a map (correspondence) between a pair of meshes from the FAUST dataset,
and create an image that visualizes this correspondence.

This code was written by Etienne Corman and modified by Maks Ovsjanikov.

Note that there are several differences between this code and the one used
to produce the results reported in the paper:

1) This implementation uses a quasi-Newton optimization strategy, rather
than writing down the linear system Ax = b explicitly. This significantly
simplifies the code, but the solution relies on an iterative solver (minConf).

2) Rather than automatically detecting segment correspondences, this code
uses a small set of pointwise landmarks. Replacing them with segment-wise
correspondences can be done by modifying lines 44 and 52 in test_compute_faust.m.

3) Neither the exact choice of descriptors, nor their parameters (values of WKS)
is the same as used in the paper. The parameters provided in this code should 
give reasonable results but are in no way guaranteed to be optimal.

4) The map evaluation criteria (final output of the code) are based on Euclidean
distances, rather than geodesic ones, and only give an aggregate mean error,
rather than complete error plots.


