function [ LacImRNA, TetRmRNA, LacI, TetR, varargout ] = generate_data_explicit_inputs(p,Y0,tspan,inputs,varargin)
%% inputs is defined by 3 column vectors, each tuple being time, ATC and IPTG
% the syntax is like in matlab function stairs: at that time, ATC and IPTG
% to their given values until a change or the end
% tspan is a vector of time points that can be identical to or different
% from inputs(:,3)
% ODE and SSA simulations are supported
% tspan is expected in seconds, unlike generate_data that expects time in
% minutes

%% Parsing inputs
ip = inputParser;
addRequired(ip,'p',@isstruct);
addRequired(ip,'Y0',@isvector);
addRequired(ip,'tspan',@isvector);
addRequired(ip,'inputs',@ismatrix);
addParameter(ip,'simulation_method','ODE',@ischar);
parse(ip,p,Y0,tspan,inputs,varargin{:});

%detects if tspan shorter than last planned change in media (ie truncated
%experiment) and if yes adapt inputs vector
if inputs(end,1)>tspan(end)
    disp('warning: experiment shorter than planned');%should not happen. can be erased
    keyboard
    ilast_change= find(inputs(:,1)>tspan(end),1,'first');
    %truncate input vector until moment of last effective change
    inputs= inputs(1:ilast_change-1,:); 
end

%identify indices at which inputs change (from c-1 to c).
c_atc= find(diff(inputs(:,2))); 
c_iptg= find(diff(inputs(:,3))); 
c_any= union(c_atc, c_iptg); % is sorted
c= [1; c_any+1]; %index of time of change

switch ip.Results.simulation_method
    case 'ODE'
        opts = odeset('NonNegative', 1:5);
    case 'SSA'
        [stoich_matrix, propensities, reactions, species_names] = toggle_props_stoich();
        iptgbefore = Y0(end-1); % IPTG is always last
        atcbefore = Y0(end); % IPTG is always last
end

Y = [];

for ind = 1:length(c) %c contains at least 1
    %simulates from c(ind) to c(ind+1) or end
    c(ind);
    ispan_beg = find(tspan >= inputs(c(ind),1),1,'first'); %should also work if tspan(1)>inputs(c(ind),1)
    if ind < length(c) %not the last change
        t_end= inputs(c(ind+1),1);
        ispan_term= find(tspan < t_end,1,'last'); %should also work if tspan(end)<inputs(c(ind+1),1)
    else
        t_end= tspan(end);
        ispan_term= length(tspan);
    end
    tspan_c= unique([inputs(c(ind),1) tspan(ispan_beg:ispan_term) t_end]);
    switch ip.Results.simulation_method
        case 'ODE'
            inputs(c(ind),2:3);
            [~,Y_short] = ode15s(@toggle_derivative_sim,tspan_c/60,Y0,opts,inputs(c(ind),2:3)', p);
            if numel(tspan_c) == 2 % With a tspan of only two elements. ode15s spits out the intermediary tsteps so here i get rid of them
                Y_short = [Y_short(1,:); Y_short(end,:)];
            end
        case 'SSA'
            % Precalculate internal IPTG concentration
            % assumes time is in minutes 
            [time_iptg, iptgdelayed]= iptgdelayfcn(inputs(c(ind),3),tspan_c/60,p,iptgbefore);
            [time_atc, atcdelayed]= atcdelayfcn(inputs(c(ind),2),tspan_c/60,p,atcbefore);
            rate_params.pre_comp_iptg_v= iptgdelayed;
            rate_params.pre_comp_iptg_t= time_iptg;
            rate_params.pre_comp_atc_v= atcdelayed;
            rate_params.pre_comp_atc_t= time_atc;
            %rate_params.atc= inputs(c(ind),2);
            rate_params.reactions= reactions;
            rate_params.species_names= species_names;
            rate_params.p= p;
                        % assumes time is in minutes; 
            % MAX_OUTPUT_LENGTH= 1e6 by default; too short
            [~, Y_short]= firstReactionMethod(stoich_matrix,propensities, ...
                tspan_c/60, round(Y0(1:(end-2))),rate_params,[],1e7);
            % Add IPTG to the result
            %interp1 takes a lot of time. replace by more efficient code
            Y_short(:,end+1) = interp1(time_iptg,iptgdelayed,tspan_c/60,'linear','extrap');            
            Y_short(:,end+1) = interp1(time_atc,atcdelayed,tspan_c/60,'linear','extrap');
            iptgbefore= iptgdelayed(end);
            atcbefore= atcdelayed(end);
    end
    
    % re-init Y0
    Y0= Y_short(end,:);
        
    icat_beg= 1;
    % values at instants of change should not appear, unless they are part of tspan  
    if tspan(ispan_beg)~=inputs(c(ind),1)
        icat_beg= 2; 
    end
    icat_term= size(Y_short,1)-1;
    if ispan_term == length(tspan)
        icat_term= size(Y_short,1);
    end
    % Concatenate:
    Y = cat(1,Y,Y_short(icat_beg:icat_term,:));
end
LacImRNA= Y(:,1);
TetRmRNA= Y(:,2);
LacI= Y(:,3);
TetR= Y(:,4);
if nargout == 6 % If user asks for calculated internal iptg concentration
    varargout{1} = Y(:,5);
    varargout{2} = Y(:,6);
end
end

