function handles = draw_nullclines(inputs,p,LacI_vector,TetR_vector)
% This function simply plots the nucclines of the toggle switch according
% to the current value of the iptg & atc input. It uses the get_nullclines
% function
    hold on
    [LacI_vector, TetR_vector, nLacI_vector, nTetR_vector] = get_nullclines(inputs, p, LacI_vector, TetR_vector);
    % Display it:
    handles.TetRNlcln = plot(nLacI_vector, TetR_vector, 'Color', [1 .2 .2], 'LineWidth', 3);
    handles.LacINlcln = plot(LacI_vector, nTetR_vector, 'Color', [.2 1 .2], 'LineWidth', 3); 
end

