function f = dijkstra(M, p)
%f = DIJKSTRA(M, p) 
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
        lengths = vecnorm(M.VERT(neigs, :) - M.VERT(u, :), 2, 2);
        % Update the distance, if needed
        alt_dist = dist(u) + lengths;
        I = find(alt_dist < dist(neigs));
        dist(neigs(I)) = alt_dist(I);
        prev(neigs(I)) = u;
%         for i = 1:length(neigs)
%             if alt_dist(i) < dist(neigs(i))
%                 dist(neigs(i)) = alt_dist(i);
%                 prev(neigs(i)) = u;
%             end
%         end
    end
    
    f = dist;
end

