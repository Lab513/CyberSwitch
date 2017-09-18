%% Optimize parameters values by fitting them to the experimental calibration data
%starting with p_core_parameters

clear all
%close all
%% Add the core functions to the current path:
addpath('..');
addpath(genpath(['..' filesep 'ModelSimulation']));
addpath(genpath(['..' filesep 'Utilities']))
addpath(genpath(['..' filesep 'Data']))

load('parameters2');

data_to_compare = {'Calibration_1','Calibration_2','Calibration_3',...
                 'Calibration_4','Calibration_5','Calibration_6'};

%% Set limits for parameters search
limits.miptg = [2 4];
limits.matc = [2 4];
limits.nt = [2 4];
limits.nl = [2 4];


%% Define set of parameters that won't be searched for...

keep_steady= {'delta_mrnal','delta_mrnat','deltal','deltat','cil','cit','cl','ct','nl','nt','k_l','k_t','matc','miptg','crl','crt'};

Prots_minmax= [];
for ind=1:length(keep_steady)
    Prots_minmax(ind,1:2)= [p.(keep_steady{ind}) NaN];
end

for ind=1:length(keep_steady)
    p_equations.(keep_steady{ind}) = @(P_mm, q) P_mm(ind,1); 
end

if isempty(keep_steady)
    p_equations= struct([]);
end

%% Load up the datasets:

for data_set=1:length(data_to_compare)
    data_comp= data_to_compare{data_set};
    disp(['Loading up data from: ' data_comp]);
    data(data_set).name = data_comp;
    data(data_set).val = load(['../Data/Experimental/', data_comp,'/FittingData.mat']);
    data(data_set).ind = load(['../Data/Experimental/', data_comp,'/FittingInputs.mat']);
    
    % Adapt the levels with the lamp_factor
    data(data_set).val.rfpMothers=data(data_set).val.rfpMothers*data(data_set).val.lamp_mult;
    data(data_set).val.gfpMothers=data(data_set).val.gfpMothers*data(data_set).val.lamp_mult;
    % If time data starts negative, just shift it.
    data(data_set).val.timerfp = data(data_set).val.timerfp - data(data_set).val.timerfp(1);
    data(data_set).val.timegfp = data(data_set).val.timegfp - data(data_set).val.timerfp(1);
end

%% Experimental fitting weights:
xpnames = data_to_compare;
 
weights = [1,1,1,1,1,1];

for data_set = 1:numel(data_to_compare)
    checksout= find(strcmp(data_to_compare{data_set},xpnames));
    if checksout
        data(data_set).val.weight = weights(checksout);
    else
        data(data_set).val.weight = 1;
    end
    disp(['XP ' data_to_compare{data_set} ' - weight = ' num2str(data(data_set).val.weight)])
end


figure()
plot_mean_optimised_results(p, data, 'RunSSA',true);

keyboard
for n=1:1
%% CMAES
%Parameters appearing in p_equations won't be searched for
[p_std,p_str]= params2vector(p, fieldnames(p_equations));
p_search= ones(numel(p_std),1);
%LBounds and UBounds are defined relatively to the reference parameter
%values
opts = cmaes();
opts.MaxFunEvals= 10;%10e3;
opts.Restarts = 0;

opts.PopSize= (4+floor(3*log(length(p_std))))*3; %3 times standard population size
[opts.LBounds, opts.UBounds] = generatelimits(p,limits, fieldnames(p_equations));

% Generate different filenames for each task in the parfor loop:
numParallel = 6;
maxWorkers = 6;
%Greg
numParallel = 1; 
maxWorkers = 1;

for l = 1:numParallel
    optsCell{l} = opts;
    optsCell{l}.SaveFilename = ['variablescmaes_par' num2str(l) '.mat'];
end
 
totaltime = tic();
disp([datestr(now) ' CMA-ES loop starting'])

for (l=1:numParallel)  
    disp(['Starting parallel cmaes #' num2str(l) ', intermediate results will be saved in: ' optsCell{l}.SaveFilename])
    init_cost= mycost(p_search, [],p_std, p_str, data, p_equations, Prots_minmax)
    disp(['Starting cost: ' num2str(init_cost)]);
    %p_min= p_search; cmin= 0;
    [p_min, cmin] = cmaes('mycost', p_search, p_search/3, optsCell{l},[], p_std, p_str, data, p_equations, Prots_minmax);
    p_opt(1,:)= p_min;
    c_min(1,:)= cmin;
end


disp([datestr(now) ' CMA-ES loop over'])
toc(totaltime);

%% Load back from files if necessary:
if ~exist('p_opt','var') % This usually means the fitting was stopped before it converged, which means we have to load it back up from the files
    toload = 1:numParallel;

    for ind1 = 1:numel(toload)
        disp(['Loading from file: ' optsCell{toload(ind1)}.SaveFilename]);
        loaded = load(optsCell{toload(ind1)}.SaveFilename);
        p_opt(ind1,:) = loaded.bestever.x;
        c_min(ind1,:) = loaded.bestever.f;
    end
end
        
%% Convert back one of them and display it: (indFit to choose which one)
indFit = 1;
p_optimized(n,:) = vector2params(p_opt(indFit,:).*p_std,p_str, p_equations, Prots_minmax);

%% plot the results of the optimisation (See also AllDisplayScripts.m)
figure(n)
cla
plot_mean_optimised_results(p_optimized(n,:), data);


%% Plot cost at every eveluation
load 'outcmaesfit.dat';
cst_p_optimized(:,n)=outcmaesfit(:,6); % cost
eval_p_optimized(:,n)=outcmaesfit(:,2); %evaluation
figure(10+n)
plot(eval_p_optimized(:,n),cst_p_optimized(:,n))
xlabel('Evaluation')
ylabel('Cost')

end

%% Optimized parameters saved as "p_optimizedd" in "CMAES_output_p_optimized"
%save('CMAES_output_p_optimized','p_optimized','cst_p_optimized','eval_p_optimized')
 
