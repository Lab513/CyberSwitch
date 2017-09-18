function [ LacImRNA, TetRmRNA, LacI, TetR, iptg, atc ] = generate_data(p,Y0,tspan,inputs,varargin)
%generate_data - Main toggle switch simulation function, for ODEs and SSA
%This function simulates the toggle switch models described by either
%toggle_derivative_sim (ODE) or toggle_props_stoich (SSA) with the
%parameters p found through the fitting over the timepoints given in the
%variable tspan. The inputs can be given either under the form of an
%"event matrix" or sampled value matrix.

% [ LacImRNA, TetRmRNA, LacI, TetR ] = generate_data(p,Y0,tspan,inputs)
% returns the simulated values of all 4 species depending on the
% parameters structure p, the initial state of the species (including IPTG) 
% Y0, and the inputs as an Nx2 array of N inputs values sampled every
% minute (ATC first, IPTG second)
% 
% [ LacImRNA, TetRmRNA, LacI, TetR , iptg] = generate_data(p,Y0,tspan,inputs)
% Additionaly returns the IPTG values returned by the simulation over tspan
%
% [ ... ] = generate_data(p,Y0,tspan,inputs,'inputs_are_events',inputs_format_str)
% specifies whether the format of the inputs is events or not
% (inputs_format_str = 'yes' or 'no'). Events are described in a Nx3 matrix
% with the atc and iptg states in the first two columns and then the time
% point at which each state stops of each state on the third column.
% default format is sampled inputs value every minute in an Nx2 matrix.
% 
% [ ... ] = generate_data(p,Y0,tspan,inputs,'simulation_method',method_str)
% specifies whether to run deterministic ODE simulations or stochasitc SSA
% simulations. (method_str = 'ODE' or 'SSA')


    %% Parsing inputs
    ip = inputParser;
    addRequired(ip,'p',@isstruct);
    addRequired(ip,'Y0',@isvector);
    addRequired(ip,'tspan',@isvector);
    addRequired(ip,'inputs',@ismatrix);
    addParameter(ip,'inputs_are_events','no',@ischar);
    addParameter(ip,'simulation_method','ODE',@ischar);
    addParameter(ip,'verbosity',0,@isnumeric);
    parse(ip,p,Y0,tspan,inputs,varargin{:});
    
    % If necessary reconstruct events variable:
    switch ip.Results.inputs_are_events
        case 'no'
            events = reconstructevents(inputs);
        case 'yes'
            events = inputs;
    end
    
    switch ip.Results.simulation_method
        case 'ODE'
            % ODE options
            %Greg: commented tolerance changes: keep the default
            %opts = odeset('RelTol',1e-1,'AbsTol',1e-2);
            opts = odeset('NonNegative', 1:5);
        case 'SSA'
            % Generate stochastic model:
            [stoich_matrix, propensities, reactions, species_names] = toggle_props_stoich('verbosity',ip.Results.verbosity);
            iptgbefore = Y0(5); %
            atcbefore = Y0(6);
    end
    
    %% Simulating:
    Y = [];
    for ind1 = 1:size(events,1)
        % Recompute time-spans
        tspan_beg = find(tspan >= events(ind1,3),1,'first');
        if isempty(tspan_beg)
            break;
        end
        if ind1 < size(events,1) && events(ind1+1,3) < tspan(end)
            tspan_end = find(tspan < events(ind1+1,3),1,'last');
            last_element = events(ind1+1,3); % Will be used for Y0
        else
            tspan_end = numel(tspan);
            last_element = tspan(end)+1; % Just a small time horizon to run the simulation for the last event
        end
        tspan_short = [tspan(tspan_beg:tspan_end) last_element]-events(ind1,3);
        % Run the simulation:
        
        switch ip.Results.simulation_method
            case 'ODE'
                if ~isfield(p,'maturationt')
                    funname= @toggle_derivative_sim;
                else
                    funname= @toggle_derivative_sim_maturation;
                end
                
                [~,Y_short] = ode15s(funname,tspan_short ,Y0,opts,events(ind1,:)', p);
                if numel(tspan_short) == 2 % With a tspan of only two elements. ode15s spits out the intermediary tsteps so here i get rid of them
                    Y_short = [Y_short(1,:); Y_short(end,:)];
                end
            case 'SSA'
                % Precalculate internal iptg concentration
                [time_iptg, iptgdelayed] = iptgdelayfcn(events(ind1,2),tspan_short, p, iptgbefore);
                [time_atc, atcdelayed] = atcdelayfcn(events(ind1,1),tspan_short, p, atcbefore);
                rate_params.pre_comp_iptg_v = iptgdelayed;
                rate_params.pre_comp_iptg_t = time_iptg;
                rate_params.pre_comp_atc_v = atcdelayed;
                rate_params.pre_comp_atc_t = time_atc;
                %rate_params.atc = events(ind1,1);
                rate_params.reactions = reactions;
                rate_params.species_names = species_names;
                rate_params.p = p;               
                % Run Gillespie's SSA
                [~, Y_short] = firstReactionMethod( stoich_matrix, propensities, tspan_short, round(Y0(1:(end-2))),...
                                         rate_params);                                     

                              
                % Add calculated IPTG to the result
                %warning('IPTG is not last position if maturation is used');
                Y_short(:,5) = interp1(time_iptg,iptgdelayed,tspan_short,'linear','extrap');
                iptgbefore = iptgdelayed(end);
                % Add calculated aTC to the result
                %warning('aTC is not last position if maturation is used');
                Y_short(:,6) = interp1(time_atc,atcdelayed,tspan_short,'linear','extrap');
                atcbefore = atcdelayed(end);
        end
        
        % Concatenate:
        Y = cat(1,Y,Y_short(1:(end-1),:));   
        % re-init Y0
            Y0 = Y_short(end,:);
    end

    
    LacImRNA=Y(:,1);
    TetRmRNA=Y(:,2);
    LacI=Y(:,3); %in case of model with maturation, this is not LacI but rather matured RFP
    TetR=Y(:,4); %in case of model with maturation, this is not TetR but rather matured GFP
    iptg=Y(:,5);
    atc=Y(:,6);
end

function events = reconstructevents(inputs)
% This function recreates an event-based inputs variable that is more
% adapted to the new ODE computation we implemented

[~,c,~] = find(diff(inputs,1,2));
c = unique([0; c]);
events = [inputs(:,c+1)' c];

end
