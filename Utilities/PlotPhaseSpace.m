function [handles] = PlotPhaseSpace(input, p, varargin)

% Inputs parsing:
    ip = inputParser;

    addRequired(ip,'input',@isnumeric);
    addRequired(ip,'p',@isstruct);
    addParameter(ip,'LacImax',4e3,@isscalar);
    addParameter(ip,'TetRmax',2.5e3,@isscalar);
    addParameter(ip,'LacImin',0,@isscalar);
    addParameter(ip,'TetRmin',0,@isscalar);
    addParameter(ip,'bkgdDiv',10,@isscalar);
    addParameter(ip,'Nullclines','on',@(x) any(strcmp(x,{'on','off'})));
    addParameter(ip,'PhsPort_step',18,@isscalar);
    addParameter(ip,'NullCl_step',300,@isscalar);
    % Phase Portrait:
    ip.addParameter('BKGD_magnification_factor',3,@isscalar);
    ip.addParameter('BKGD_image','intensity',@(x) any(strcmp(x,{'intensity','orientation','rfp-gfp','none'})))
    ip.addParameter('quiver','on',@(x) any(strcmp(x,{'on','off'})))
    ip.addParameter('lines','off',@(x) any(strcmp(x,{'on','off'})))
    
    parse(ip,input,p,varargin{:});
    LacImax = ip.Results.LacImax;
    TetRmax = ip.Results.TetRmax;
    LacImin = ip.Results.LacImin;
    TetRmin = ip.Results.TetRmin;
    bkgdDiv = ip.Results.bkgdDiv;
    plotNullcl = ip.Results.Nullclines;
    PhsPort_step = ip.Results.PhsPort_step;
    NullCl_step = ip.Results.NullCl_step;

%% Plot Nullclines and Phase-space

% bkgd: (This technique with imref2d lowers memory problems)
% bkgd = getBackgroundRG(LacImax/bkgdDiv, TetRmax/bkgdDiv);
% RI = imref2d(size(bkgd));
% RI.XWorldLimits = [0 LacImax];
% RI.YWorldLimits = [0 TetRmax];
% handles.bkgd = imshow(bkgd,RI);
set(gca,'YDir','normal')
hold on
FontSize= 18;
% Phase-portrait:
LacI_vector= linspace(LacImin,LacImax,PhsPort_step);
TetR_vector= linspace(TetRmin,TetRmax,PhsPort_step);
draw_phase_space(input',p, LacI_vector,TetR_vector, ...
    'BKGD_magnification_factor',ip.Results.BKGD_magnification_factor, ...
    'BKGD_image',ip.Results.BKGD_image, ...
    'quiver',ip.Results.quiver, ...
    'lines',ip.Results.lines ...
    );

% Nullclines:
if strcmp(plotNullcl,'on')
    LacI_vector= linspace(LacImin,LacImax,NullCl_step);
    TetR_vector= linspace(TetRmin,TetRmax,NullCl_step);
    handles.nullclines = draw_nullclines(input',p, LacI_vector,TetR_vector);

    %Plots nullcline intersections
    [X0,Y0]= get_intersections(input,p, LacI_vector,TetR_vector);
    if numel(X0)== 1
        plot(X0,Y0,'o','MarkerFaceColor',[0.5,0.5,0.5],'MarkerSize',8)    
    elseif numel(X0)== 3
        u= find(X0<max(X0) & Y0<max(Y0));
        plot(X0(u),Y0(u),'o','MarkerFaceColor',[0,0,0],'MarkerSize',8)
        s= setdiff(1:3,u);
        plot(X0(s),Y0(s),'o','MarkerFaceColor',[0.5,0.5,0.5],'MarkerSize',8)
    else
    end
end

% Re-set view
xlim([LacImin LacImax]);
ylim([TetRmin TetRmax]);
xlabel('LacI-RFP', 'Color','r','FontSize',FontSize, 'FontWeight','Bold');
ylabel('TetR-GFP', 'Color','g','FontSize',FontSize, 'FontWeight','Bold');
set(gca,'LineWidth',2, 'FontSize',FontSize);%set(ax(1),'Ycolor',[1 0 0]);
set(gca,'LineWidth',2, 'FontSize',FontSize);%set(ax(2),'Ycolor',[0 1 0]);
