function ftri = func_on_tris(M, f)
    ftri = sum(f(M.TRIV), 2) ./ 3;
end

