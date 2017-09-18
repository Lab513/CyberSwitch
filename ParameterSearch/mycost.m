function [total_cost] = mycost(p_search, ~,p_std, p_str, data, p_eq, P_mm)
%iteration on all the data.m files in data_for_optimisation
    % Goal of this part is to avoid getting negative parameters
            
    % Reformat the parameters:
    p = vector2params(p_search.*p_std', p_str, p_eq, P_mm);
    
    % initialisation of the cost
    total_cost=0;
    
    % Run the loop on all data:
    for ds=1:numel(data)
        
        % Compute starting point based on overnight culture:
        %Y0=startingpoint_value( p ,data(ds).ind.inputsbefore , data(ds).val);
        Y0=startingpoint_ode( p ,data(ds).ind.inputsbefore, data(ds).val);
        
        % Run the simulation with the new paramters:
        [~,~,LacI,TetR] = generate_data(p,Y0,data(ds).val.timerfp/60,data(ds).ind.inputs); 

        % Compute cost:
            total_cost = total_cost + cost_fcn( data(ds).val.rfpMothers, LacI, data(ds).val.timerfp(end)/60, data(ds).val.weight);
            total_cost = total_cost + cost_fcn( data(ds).val.gfpMothers, TetR, data(ds).val.timerfp(end)/60, data(ds).val.weight);
    end
end