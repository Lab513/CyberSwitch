function [timepoints, iptgdelayed] = iptgdelayfcn(i_ext,tspan, p, iptgbefore)
%time is assumed to be in minutes

% iptg1storder = @(t,y,i_ext,p) (i_ext - y)/(p.iptgdelay1/60);
iptg1storder = @(t,y,i_ext,p) max(((i_ext-y)/(p.iptgdelay1/60)),0)-max(((y-i_ext)/(p.iptgdelay2/60)),0);
opts = odeset('NonNegative',1);
[timepoints, iptgdelayed] = ode15s(iptg1storder,tspan,iptgbefore,opts,i_ext,p);



        