function [p] = vector2params(p_vec, p_str, p_eq, P_mm)
% A function to convert parameters structure to a vector once and for all

for ind1 = 1:numel(p_str)
    p.(p_str{ind1}) = p_vec(ind1);
end

% Processing equations:
p_eqstr = fieldnames(p_eq);
for ind1 = 1:numel(p_eqstr)
    p.(p_eqstr{ind1}) = p_eq.(p_eqstr{ind1})(P_mm, p);
end