function h = hacksurf(x,y,t,LineWidth)

x = x';
y = y';
% col = 1:numel(x);

z = zeros(size(x))*max(x);
h = surface([x;x],[y;y],[z;z],[t;t],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',LineWidth);