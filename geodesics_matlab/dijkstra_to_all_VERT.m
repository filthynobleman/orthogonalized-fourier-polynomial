function D = dijkstra_to_all_VERT(shape, sources)
    S = shape;
    
    if(size(sources,2) > size(sources,1))
        sources = sources';
    end
    
    if(size(sources,2) > 1)
        error('sources must be stored in a vector');
    end
    
    D = comp_geodesics_to_all(double(S.VERT(:,1)), double(S.VERT(:,2)), double(S.VERT(:,3)), ...
                              double(S.TRIV'), sources, 1);
end