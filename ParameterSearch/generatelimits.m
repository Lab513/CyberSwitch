function [LBounds, UBounds] = generatelimits(p, lims, donotprocess, varargin)
% [LBounds, UBounds] = generatelimits(p, lims, donotprocess, varargin)
%
% This function generates the lower and upper bounds variables LBounds and UBounds that will be
% used by CMA-ES based on the parameters structure variable p and the limis defined
% in the structure variable lims. The parameters that have no limits defined in lims
% will be limited to default [1e-3, 1e3]. If some limits should not be processed, the
% donotprocess cell containing their names should be provided, otherwise
% an empty cell can replace it. One optional parameter can be provided
% to change the default limits (must be a two elements vector like [1e-3
% 1e3])

%this default not used anymore
% if nargin >= 4
%     default = varargin{1};
% else
%     default = [1e-3 1e3];
% end

p_str = setdiff(fieldnames(p),donotprocess);
l_str = fieldnames(lims);
for ind1 = 1:numel(p_str)
    if any(strcmp(p_str{ind1}, l_str))
        LBounds(ind1,1) = lims.(p_str{ind1})(1)/p.(p_str{ind1});
        UBounds(ind1,1) = lims.(p_str{ind1})(2)/p.(p_str{ind1});
    else
        LBounds(ind1,1) = 1/5;% default(1); % All parameters must be > 0. I do not use 0 but a very low value to avoid NaN problems (div by 0 etc...)
        UBounds(ind1,1) = 5; % default(2); % Same for upper bounds...
    end
end
