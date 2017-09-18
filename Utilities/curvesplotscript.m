function curvesplotscript(XP,varargin)
%% Input processing:

possiblemethods = {'levels','ratio','levelslog','color2D'};
possiblewhich = {'all','controlled','mean'};
ip = inputParser();
addRequired(ip,'XP',@isstruct);
addParameter(ip,'method','levels',@(x) any(strcmp(x,possiblemethods)));
addParameter(ip,'which','all',@(x) any(strcmp(x,possiblewhich)) || isstrprop(x,'digit'));
addParameter(ip,'LacImax',3.1e3,@isscalar)
addParameter(ip,'TetRmax',1.5e3,@isscalar)
parse(ip,XP,varargin{:});

%if ~any(strcmp(fieldnames(XP.gfp_ctrl.controller),'objective'))%undefined for BangBang controllers
if ~isempty(XP.type) && strcmp(XP.type,'control') && isempty(XP.gfp_ctrl.controller)
    XP.gfp_ctrl.controller.objective= [350,Inf];
    XP.rfp_ctrl.controller.objective= [750,Inf];
end

LacImax = ip.Results.LacImax;
TetRmax = ip.Results.TetRmax;

%% Plotting:
if ~strcmp(ip.Results.method,'color2D')
    switch ip.Results.method
        case {'levels', 'levelslog'}
            if strcmp(XP.type,'control')
                c2c = XP.rfp_ctrl.celltocontrol;
                objR= XP.rfp_ctrl.controller.objective(1)*ones(length(XP.timepoints),1);
                objG= XP.gfp_ctrl.controller.objective(1)*ones(length(XP.timepoints),1);
            else
                objR= []; objG=[];
            end
            switch ip.Results.which
                case 'all'
                    [ax,hr,hg] = plotyy(XP.timepoints/3600,[XP.rfp objR],XP.timepoints/3600,[XP.gfp objG]);
                    set(hr,'Color',[1 0 0]); %pure red
                    set(hg,'Color',[0 1 0]); %pure green
                    if strcmp(XP.type,'control')
                        %set(hr(setdiff(1:numel(hrf),c2c)),'Color',[1 .6 .6])
                        set(hr([c2c end]),'LineWidth',2); %darker blueish red
                        set(hr(end),'LineStyle','--'); 
                        set(hg([c2c end]),'LineWidth',2); %darker blueish green
                        set(hg(end),'LineStyle','--');
                    end

                case 'mean'
                    [ax,hr,hg] = plotyy(XP.timepoints/3600,[mean(XP.rfp,2) objR],XP.timepoints/3600,[mean(XP.gfp,2) objG]);
                    set(hr,'LineWidth',2,'Color',[1 0 0]);
                    set(hg,'LineWidth',2,'Color',[0 1 0]);
                case 'controlled'
                    [ax,hr,hg] = plotyy(XP.timepoints/3600,[XP.rfp(:,c2c) objR],XP.timepoints/3600,[XP.gfp(:,c2c) objG]);
                    set(hr,'LineWidth',2,'Color',[1 0 0]);
                    set(hg,'LineWidth',2,'Color',[0 1 0]);
                otherwise
                    tp = str2num(ip.Results.which);
                    [ax,hr,hg] = plotyy(XP.timepoints/3600,[XP.rfp(:,tp) objR],XP.timepoints/3600,[XP.gfp(:,tp) objG]);
                    set(hr,'LineWidth',2,'Color',[1 0 0]);
                    set(hg,'LineWidth',2,'Color',[0 1 0]);
            end

            ylabel(ax(1),'LacI-RFP')
            ylabel(ax(2),'TetR-GFP')
            ylim(ax(1),[0 LacImax]);
            ylim(ax(2),[0 TetRmax]);
            xlim(ax(1),[0 XP.timepoints(end)/3600])
            xlim(ax(2),[0 XP.timepoints(end)/3600])     
            set(ax(1),'YTick',0:500:LacImax);set(ax(1),'Ycolor',[1 0 0]);
            set(ax(2),'YTick',0:500:TetRmax);set(ax(2),'Ycolor',[0 1 0]);
            set(ax(1),'FontSize',12);
            set(ax(2),'FontSize',12);

            if strcmp(ip.Results.method,'levelslog')
                set(ax(1),'Yscale','log');
                set(ax(2),'Yscale','log');
            end

        case 'ratio'
            %title('Fluorescence ratio')
            % RFP/GFP
            if strcmp(XP.type,'control') 
                c2c = XP.rfp_ctrl.celltocontrol;
                objctv = XP.rfp_ctrl.controller.objective(1)/XP.gfp_ctrl.controller.objective(1);
            end
            switch ip.Results.which
                case 'all'
                    hrf = plot(XP.timepoints/3600,XP.rfp./XP.gfp,'Color',[.85 .33 .1]);
                case 'controlled'
                    hrf = plot(XP.timepoints/3600,XP.rfp(:,c2c)./XP.gfp(:,c2c),'Color',[.85 .33 .1]);
                case 'mean'
                    hrf = plot(XP.timepoints/3600,mean(XP.rfp,2)./mean(XP.gfp,2),'Color',[.85 .33 .1]);
                otherwise
                    tp = str2num(ip.Results.which);
                    hrf = plot(XP.timepoints/3600,XP.rfp(:,tp)./XP.gfp(:,tp),'Color',[.85 .33 .1]);
            end
            if strcmp(XP.type,'control') 
                set(hrf(c2c),'LineWidth',2)
                set(hrf(setdiff(1:numel(hrf),c2c)),'Color',[.93 .69 .13])
                line(xlim,[objctv objctv ],'Color','k')
            end
            set(gca,'Yscale','log');
            ylabel('LacI-RFP / TetR-GFP')
            ylim([1e-2 1e2]);  
            xlim([0 XP.timepoints(end)/3600])
            set(gca,'FontSize',12);
    end
    xlabel('Time (h)')
    set(gca,'XTick',0:5:24);
else
    % Smooth it first:
    span = 20/size(XP.rfp,1);
    for ind1 = 1:size(XP.rfp,2)
        x(:,ind1) = smooth(XP.rfp(:,ind1),span,'sgolay',3);
        y(:,ind1) = smooth(XP.gfp(:,ind1),span,'sgolay',3);
    end
    if strcmp(XP.type,'control') 
        c2c = XP.rfp_ctrl.celltocontrol;
        objR= XP.rfp_ctrl.controller.objective(1)*ones(length(XP.timepoints),1);
        objG= XP.gfp_ctrl.controller.objective(1)*ones(length(XP.timepoints),1);
    end
    % Plot:
    switch ip.Results.which
        case 'all'
            error('You cannot use the 2D color plot with all cells, it''s way too heavy');
        case 'controlled'
            x = x(:,c2c);
            y = y(:,c2c);
        case 'mean'
            x = mean(x,2);
            y = mean(y,2);
        otherwise
            tp = str2num(ip.Results.which);
            x = x(:,tp);
            y = y(:,tp);
            
    end
    hacksurf(x,y,XP.timepoints/3600,4);
    xlim([0 LacImax]);
    ylim([0 TetRmax]);
    if strcmp(XP.type,'control') 
        line([objR objR],[0 1500],'Color','r');
        line([0 3100],[objG objG],'Color','g');
    end
    xlabel('LacI-RFP');
    ylabel('TetR-GFP');
    colormap('jet')
    hc = colorbar();
    ylabel(hc, 'time (h)')
end
set(gca,'Box','off');


