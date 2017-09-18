function [Y] = startingpoint_ode( p ,inputsbefore, values )
%STARTING POINT compute the concentrations of mrna and proteins after a
%night for =/= inducer concentrations. the format of the inducer is
%[ATC;IPTG] it should contain only 1 value for each inducer.

    opts = odeset('NonNegative',1:5);
    %Greg: changed to avoid issues with bistable systems (the wrong
    %configuration can be found if inputsbefore = [0,0]
    %y0= [0;0;0;0;0];
    
    Y0= startingpoint_value(p, inputsbefore, values);
    funname= @toggle_derivative_sim;
    [~,Y] = ode15s(funname,[0 4800],Y0,opts,inputsbefore',p);
    Y=Y(end,:);
    
end