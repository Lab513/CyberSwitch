function [p_vec, p_str] = params2vector(p, donotprocess)
% A function to convert parameters structure to a vector once and for all

p_str = setdiff(fieldnames(p), donotprocess);

for ind1 = 1:numel(p_str)
    p_vec(ind1) = p.(p_str{ind1});
end