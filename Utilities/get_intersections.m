function [ X0,Y0 ] = get_intersections(inputs,p, LacI_vector,TetR_vector)
    [ LacI_vector, TetR_vector, nLacI_vector, nTetR_vector]=get_nullclines( inputs,p, LacI_vector,TetR_vector);
    [X0,Y0] = intersections(LacI_vector, nTetR_vector, nLacI_vector, TetR_vector);
end
