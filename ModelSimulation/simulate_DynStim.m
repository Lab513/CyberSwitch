function XP = simulate_DynStim(p,stim_iptg,stim_atc,time_iptg,time_atc,repeats,release,time_release,numcells)


% Generate stimulation inputs:
inputs_valve = [stim_iptg time_iptg; stim_atc time_atc]';
time_events = [0 cumsum(repmat(inputs_valve(3,:),1,repeats))];
events = cat(1,repmat(inputs_valve(1:2,:),1,repeats),time_events(1:(end-1)));
%Add release at the end:
events = cat(2, events, [ release time_events(end)]');
tspan = 0:5:(events(3,end)+ time_release );

% Generate data:
inputs_before= [0 1];
Y0= startingpoint_ode(p, inputs_before, struct('rfpMothers',0,'gfpMothers',0));
parfor indCell = 1:numcells
    tic
    [~,~,LacI(indCell,:),TetR(indCell,:),~,~] = generate_data(p,Y0,tspan,events,'simulation_method','SSA');
    fprintf('Simulated %d/%d cells\n',indCell,numcells);
    toc
end

events_reform = events([2 1 3],:)';
events_reform(:,3) = cumsum([diff(events_reform(:,3)); time_release]*60);

XP= struct('rfp',LacI, 'gfp',TetR, 'timepoints',tspan.*60, 'media', events_reform, 'type', 'open-loop');