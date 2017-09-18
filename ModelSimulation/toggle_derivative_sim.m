function dydt = toggle_derivative_sim(t,y,inputs,p)
%% Deterministic Model of the Toggle Switch

%% Handling inputs
% Set the inducer values depending on the format of the inputs variable:
if size(inputs,2) == 1 % We got only one set of values (usually for overnight calculations)
        atc_ext = inputs(1);
        iptg_ext = inputs(2);
elseif floor(t)+1 <= size(inputs,2) % I use floor + 1 instead of round because it avoids access to element 0 cases (although that shouldn't be too much of a pb anyways)
        atc_ext =  inputs(1,floor(t)+1);
        iptg_ext = inputs(2,floor(t)+1);
else
       warning('Toggle model: The inputs variable is not long enough for the requested simulation time. Using last inducers value in the inputs variable')
       atc_ext =  inputs(1,end);
       iptg_ext = inputs(2,end);
end




%% Choose between quiver and ODE mode:

switch size(y,3)
    case 1 % ODE / normal mode
        dydt= zeros(5,1);
        dydt(1) = p.crl + p.cil * hill_func(y(4) * hill_func(y(6),p.katc,p.matc),p.k_t,p.nt) - p.delta_mrnal*y(1);%dlaciM(y(4),atc,y(1),p);
        dydt(2) = p.crt + p.cit * hill_func(y(3) * hill_func(y(5),p.kiptg,p.miptg),p.k_l,p.nl) - p.delta_mrnat*y(2);%dtetrM(y(3),y(5),y(2),p);
        dydt(3) = p.cl*y(1) - p.deltal*y(3);
        dydt(4) = p.ct*y(2) - p.deltat*y(4);
%                 dydt(5) = (iptg_ext - y(5))/(p.iptgdelay1/60);
                 dydt(5) = max((iptg_ext - y(5))/(p.iptgdelay1/60),0)-max((y(5)-iptg_ext)/(p.iptgdelay2/60),0);

     dydt(6) = max((atc_ext - y(6))/(p.atcdelay1/60),0)-max((y(6)-atc_ext)/(p.atcdelay2/60),0);
    case 4 % Quiver mode
        dydt= zeros([size(y),4]);
        dydt(:,:,3) = p.cl*y(:,:,1) - p.deltal*y(:,:,3); %dlaciP(y(:,:,1),y(:,:,3),p); 
        dydt(:,:,4) = p.ct*y(:,:,2) - p.deltat*y(:,:,4); %dtetrP(y(:,:,2),y(:,:,4),p);
end
dydt= real(dydt);
end