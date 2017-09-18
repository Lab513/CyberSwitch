function [timepoints, atcdelayed] = atcdelayfcn(i_ext,tspan, p, atcbefore)
%time is assumed to be in minutes

atc1storder = @(t,y,i_ext,p) max(((i_ext - y)/(p.atcdelay1/60)),0)-max(((y-i_ext)/(p.atcdelay2/60)),0);
opts = odeset('NonNegative',1);
[timepoints, atcdelayed] = ode15s(atc1storder,tspan,atcbefore,opts,i_ext,p);
