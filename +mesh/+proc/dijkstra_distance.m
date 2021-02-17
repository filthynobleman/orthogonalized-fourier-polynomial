function f = dijkstra_distance(M, p)
%f = DIJKSTRA_DISTANCE(M, p) 
%   Approximates the geodesic distance function from a vertex p on the mesh
%   M using the Dijkstra distance on the connectivity of M.

    % Create a set of vertices
    Q = true(M.n, 1);
    % Create the vectors of current and previous distance of the vertices
    dist = inf(M.n, 1);
    dist(p) = 0;
    prev = zeros(M.n, 1);
    
    while any(Q)
        % Get the vertex with minimum distance and remove it from Q
        u = find(and(Q, (dist == min(dist(Q)))), 1, 'first');
        Q(u) = false;
        
        % Find the neighbours of u
        neigs = M.TRIV(any(M.TRIV == u, 2), :);
        neigs = unique(neigs(:));
        % Get the euclidean distance along the edges
        lengths = utils.repnorm(M.VERT(neigs, :) - M.VERT(u, :), 'rows');
        % Update the distance, if needed
        alt_dist = dist(u) + lengths;
%         dist(alt_dist < dist) = alt_dist(alt_dist < dist);
%         prev(alt_dist < dist) = u;
        for i = 1:length(neigs)
            if alt_dist(i) < dist(neigs(i))
                dist(neigs(i)) = alt_dist(i);
                prev(neigs(i)) = u;
            end
        end
    end
    
    f = dist;
end

