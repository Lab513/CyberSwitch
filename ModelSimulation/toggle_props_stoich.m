function [stoich_matrix, propensities, reactions, species_names] = toggle_props_stoich(varargin)

% Inputs parsing
    ip = inputParser;
    addParameter(ip,'verbosity',0,@isnumeric);
    parse(ip,varargin{:});

%% Stochastic Model of the Toggle Switch
% Declare species names:
species_names = {'LacImRNA', 'TetRmRNA', 'LacI', 'TetR'};

%% Declare the equations & stoichiometry: 
% "pTet" transcription of laci mrna
reactions.lacim_trscr.stoich = {'LacImRNA',+1};
reactions.lacim_trscr.prop = @(spcs,ipts,p) ...
            p.crl ...
            + p.cil .* hill_func(   spcs.TetR .* hill_func( ipts.atc, ...
                                                            p.katc, ...
                                                            p.matc), ...
                                    p.k_t, ...
                                    p.nt);
                                
% "pLac" transcription of tetr mrna
reactions.tetrm_trscr.stoich = {'TetRmRNA',+1};
reactions.tetrm_trscr.prop = @(spcs,ipts,p) ...
            p.crt ...
            + p.cit .* hill_func(   spcs.LacI .* hill_func( ipts.iptg, ...
                                                            p.kiptg, ...
                                                            p.miptg), ...
                                    p.k_l, ...
                                    p.nl);
                                
% Translation of LacI protein
reactions.LacIp_trsl.stoich = {'LacI',+1};
reactions.LacIp_trsl.prop = @(spcs,ipts,p) p.cl.*spcs.LacImRNA;

% Translation of TetR protein
reactions.TetRp_trsl.stoich = {'TetR',+1};
reactions.TetRp_trsl.prop = @(spcs,ipts,p) p.ct.*spcs.TetRmRNA;

% Degradation of laci mRNA
reactions.lacim_deg.stoich = {'LacImRNA',-1};
reactions.lacim_deg.prop = @(spcs,ipts,p) p.delta_mrnal.*spcs.LacImRNA;

% Degradation of tetr mRNA
reactions.tetrm_deg.stoich = {'TetRmRNA',-1};
reactions.tetrm_deg.prop = @(spcs,ipts,p) p.delta_mrnat.*spcs.TetRmRNA;

% Dilution of LacI proteins
reactions.LacIp_dil.stoich = {'LacI',-1};
reactions.LacIp_dil.prop = @(spcs,ipts,p) p.deltal.*spcs.LacI;

% Dilution of TetR proteins
reactions.TetRp_dil.stoich = {'TetR',-1};
reactions.TetRp_dil.prop = @(spcs,ipts,p) p.deltat.*spcs.TetR;

%% The output variables that will be used by the SSA function:
stoich_matrix = processStoichiometry(reactions,species_names,ip.Results.verbosity);
propensities = @propensites_computation;



function stoich_matrix = processStoichiometry(reactions,spcsn,verb)
% From the reactions defined in the main function above, create the
% stoichiometry matrix for the stochastic simulation algorithm:

% Initialize:
rctn_names = fieldnames(reactions);
visp('Generating stoichiometry matrix:',1,verb)
stoich_matrix = [];

% Loop through all reactions:
for ind1 = 1:numel(rctn_names)
    stoich_cel = reactions.(rctn_names{ind1}).stoich;
    stoich_vec = zeros(1,numel(spcsn));
    % Display:
    visp(sprintf(['\tReaction #' num2str(ind1) ' (' rctn_names{ind1} ') generates:\n']),1,verb)
    % Loop through all species involved:
    for ind2 = 1:size(stoich_cel,1)
        visp(sprintf(['\t\t' num2str(stoich_cel{ind2,2}) ' ' stoich_cel{ind2,1} ' molecule\n']),1,verb)
        stoich_vec(strcmp(stoich_cel{ind2,1},spcsn)) = stoich_cel{ind2,2}; % The vector for this specifiec reaction
    end
    
    % Display: 
    visp(sprintf('\t\tStoichiometry vector: '),1,verb);
    visp(num2str(stoich_vec),1,verb);
    visp(sprintf('\n'),1,verb);
    
    % The actual matrix:
    stoich_matrix = cat(1,stoich_matrix,stoich_vec);

end

% Display:
visp('Stoichiometry Matrix:',1,verb)
visp(num2str(stoich_matrix),1,verb)
visp(sprintf('\n'),1,verb);


function props = propensites_computation(spcs_vec,entire_params,t)
% This is the function that will be called by the gillespie simulation to
% compute propensities:

% Initialize:
reactions = entire_params.reactions;
spcs = reconstructspcs(spcs_vec,entire_params.species_names);
rctn_names = fieldnames(reactions);

% Inducer levels:
ipts.iptg = interp1(entire_params.pre_comp_iptg_t,entire_params.pre_comp_iptg_v,t,'linear','extrap');
%Greg on June 2 2017. Why this has not been done before?
ipts.atc = interp1(entire_params.pre_comp_atc_t,entire_params.pre_comp_atc_v,t,'linear','extrap');
%old code: ipts.atc = entire_params.atc;

% Compute propensities:
for ind1 = 1:numel(rctn_names)
    props(ind1,1) = real(reactions.(rctn_names{ind1}).prop(spcs,ipts,entire_params.p));
end



%% Utilities:
function spcs = reconstructspcs(spcs_vec,spcsn)
for ind1 = 1:numel(spcsn)
    spcs.(spcsn{ind1}) = spcs_vec(ind1);
end

function visp(msg,verbosity,verb)
if verbosity < 0
error(msg);
end
if verb >= verbosity
fprintf(msg);
end