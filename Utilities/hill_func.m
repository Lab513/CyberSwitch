function [ output_args ] = hill_func(x,k,n) 
%HILL_FUNC Summary of this function goes here
%   Detailed explanation goes here
output_args = 1./(1+(x./k).^n);
end

