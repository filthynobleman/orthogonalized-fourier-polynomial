function f = dijkstra_matrix(M)
%f = DIJKSTRA(M, p) 
%   Approximates the geodesic distance function from a vertex p on the mesh
%   M using the Dijkstra distance on the connectivity of M.

    dist = inf(M.n, M.n);
    dist(1:M.n+1:end) = 0;
    prev = zeros(M.n, M.n, 'uint32');
    lengths = squareform(pdist(M.VERT));
    
    fprintf("Calling C routine...\n");
    f = mesh.metrics.dijkstra_internal(M.VERT, M.TRIV - 1, dist, prev, lengths);
end

