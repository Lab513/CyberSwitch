function [inputs] = generate_integrative_inputs(p,inputs_time)

Y0 = zeros(1,6); % Doesn't matter for iptg and atc

time_events = [0 cumsum(repmat(inputs_time(3,:),1,20))];
events = cat(1,repmat(inputs_time(1:2,:),1,20),time_events(1:(end-1)));
tspan = [0:1:time_events(end)];

[~,~,~,~,iptg_del,atc_del] = generate_data(p,Y0,tspan,events);

idx = find(tspan>=time_events(end-2),1,'first');
inputs = [ atc_del(idx:end) iptg_del(idx:end) ones(size(atc_del(idx:end)))]';
