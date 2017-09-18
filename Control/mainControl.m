close all
clear all

%% Path
addpath('..');
addpath(genpath(['..' filesep 'ModelSimulation']));
addpath(genpath(['..' filesep 'Utilities']));

%%control type
c_type_string= {'weakPI', 'strongPI', 'BangBang'};
c_type= c_type_string{2};

%% Parameters
load('parameters2');

% Time variables (minutes)
tstep = 5;
duration = 24*60;
curr_tp = 0;
tpts = 0;

numcells = 16;

% Controllers:
if strcmp(c_type, 'weakPI') % Parameters as in experiments
  LacI_ctrl = PIcontrollerNew(0.0330,1.38e-4,[750, Inf],[0,50],20,120);
  TetR_ctrl = PIcontrollerNew(0.025,6.9e-4,[350, Inf],[0,50],25,120);
  objective= [750 300];
elseif strcmp(c_type, 'strongPI') % Parameters as in ED figure 7
    LacI_ctrl = PIcontrollerNew(5e-2,2e-4,[750, Inf],[0,50],20,120);
    TetR_ctrl = PIcontrollerNew(2.5e-2,6.94e-4,[350, Inf],[0,50],25,120);
      objective= [750 300];
elseif strcmp(c_type, 'BangBang') % Parameters as in experiments
    LacI_ctrl = BBcontrollerNew([750, Inf],[0,50]);
    TetR_ctrl = BBcontrollerNew([350, Inf],[0,50]);
      objective= [750 300];
else
    keyboard
end

%% Plot the true control result
load(['..' filesep 'ExperimentalData' filesep 'Control' filesep c_type filesep 'XP.mat']);
figure(1);subplot(3,1,3);cla;hold on
valvesplotscript(XP);
subplot(3,1,[1 2]);cla;hold on;
curvesplotscript(XP)

%% Launch the control loop
values.rfpMothers=0; %Could be the observed values in case one reproduces a control experiment
values.gfpMothers=0;
% Inputs:
iptg = 1;
atc = 0;
iptg_mult = 1/100;
atc_mult = 1;
Y0= startingpoint_ode( p ,[atc iptg], values);
LacImRNA = Y0(end,1);
TetRmRNA = Y0(end,2);
LacI = Y0(end,3);
TetR = Y0(end,4);
iptg_del = Y0(end,5);
atc_del = Y0(end,6);

while(curr_tp<duration)    
    % Get the decision of the two controllers based on last measurements:
    iptg(end+1) = TetR_ctrl.decide(TetR',tpts*60)*iptg_mult;
    atc(end+1) = LacI_ctrl.decide(LacI',tpts*60)*atc_mult;
    
    % Compute the result of that decision:
     [Lm,Tm,L,T,I,A] =       generate_data(...
                                        p,...
                                            [LacImRNA(end),...
                                            TetRmRNA(end),...
                                            LacI(end),...
                                            TetR(end),...
                                            iptg_del(end),...
                                            atc_del(end)],...
                                        [0 tstep],...
                                        [atc(end) iptg(end) 0],...
                                        'inputs_are_events','yes',...
                                        'Simulation_method','SSA',...
                                        'verbosity',0);         
    LacImRNA(end+1) = Lm(end);
    TetRmRNA(end+1) = Tm(end);
    LacI(end+1) = L(end);
    TetR(end+1) = T(end);
    iptg_del(end+1) = I(end);
    atc_del(end+1) = A(end);
    % Update time variables:
    curr_tp = curr_tp+tstep;
    tpts(end+1) = curr_tp;
    
    disp(curr_tp/duration)
    
      
end

%% Simulate other cells:

disp(['Simulating ' num2str(numcells) ' other cells...'])
parfor (indPar=1:numcells, 6)
    [othercells(indPar).LacImRNA,...
     othercells(indPar).TetRmRNA,...
     othercells(indPar).LacI,...
     othercells(indPar).TetR,...
     ] = ...
                          generate_data(...
                                        p,...
                                        Y0, ...
                                        tpts,...
                                        [atc(2:end)' iptg(2:end)' tpts(2:end)'],...
                                        'inputs_are_events','yes',...
                                        'Simulation_method','SSA',...
                                        'verbosity',0);
                                    
    fprintf('\t Finished cell #%d...\n',indPar)
end

%% Plot it all:
% LacI
figure(4);cla;hold on
% TetR
figure(5);cla;hold on
% LacI/TetR
figure(6);cla;hold on

% Plot other cells:
%Greg: lack of image toolbox
%colors = distinguishable_colors(numcells,{'w','k'});
 for ind2 = 1:numcells
    figure(4)
    %plot(tpts./60,othercells(ind2).LacI,'color',colors(ind2,:),'linewidth',1)
    plot(tpts./60,othercells(ind2).LacI,'linewidth',1)
    figure(5)
    %plot(tpts./60,othercells(ind2).TetR,'color',colors(ind2,:),'linewidth',1)
    plot(tpts./60,othercells(ind2).TetR,'linewidth',1)
    figure(6)
    %plot(tpts./60,othercells(ind2).LacI./othercells(ind2).TetR,'color',colors(ind2,:),'linewidth',1)
    plot(tpts./60,othercells(ind2).LacI./othercells(ind2).TetR,'linewidth',1)
 end

% LacI
figure(4)
plot(tpts/60,LacI,'k','LineWidth',2)
title('LacI-RFP levels')
ylabel('LacI-RFP')
xlabel('time (h)')
line(xlim,[objective(1) objective(1)],'Color','k')
ylim([0 4e3])

% TetR
figure(5)
plot(tpts/60,TetR,'k','LineWidth',2)
title('TetR-GFP levels')
ylabel('TetR-GFP')
xlabel('time (h)')
line(xlim,[objective(2) objective(2)],'Color','k')
ylim([0 2.5e3])

% LacI/TetR
figure(6)
plot(tpts/60,LacI./TetR,'k','LineWidth',2)
title('LacI/TetR ratio')
ylabel('LacI/TetR')
xlabel('time (h)')
line(xlim,[objective(1)./objective(2) objective(1)./objective(2)],'Color','k')
set(gca,'YScale','log')
ylim([1e-2 1e2])

%% LacI-TetR
figure(7)
cla
hacksurf(LacI',TetR',tpts/60,3)
% title('LacI/TetR ratio')
ylabel('TetR')
xlabel('LacI')
ylim([0 2e3])
xlim([0 3e3])
colormap('jet')
colorbar
line(xlim,[objective(2) objective(2)],'Color','g')
line([objective(1) objective(1)],ylim,'Color','r')
% set(gca,'YScale','log')


%% inputs
figure(8); hold on
[iptgstrX, iptgstrY] = stairs([0; tpts'/60],[iptg'; iptg(end)]);
[atcstrX, atcstrY] = stairs([0; tpts'/60],[atc'; atc(end)]+125);
h = area(iptgstrX, 100*iptgstrY, 0);
set(h(1),'FaceColor',[0 1 1]);
set(h(1),'EdgeColor',[0 1 1]);
h = area(atcstrX, atcstrY, 125);
set(h(1),'FaceColor',[1 0 1]);
set(h(1),'EdgeColor',[1 0 1]);
xlim([0 tpts(end)/60])
ylim([-5 250])

% Plot max and min lines:
line([xlim],[0     0],'Color','k')
line([xlim],[100 100],'Color','k')
line([xlim],[125 125],'Color','k')
line([xlim],[225 225],'Color','k')
% Cosmetics:
set(gca,'XTick',[0 5 10 15 20 25 30 40 50 60 ]);
set(gca,'YTick',[0 100 125 225]);
set(gca,'YTickLabel',{'0%' '100%' '0%' '100%'});
set(gca,'Box','off');
xlabel('time (h)')
ylabel('IPTG                             aTC')


%% Phase-space over time:
if 1 %used to skip the plot
    for ind1 = 1:numel(tpts)
        figure(3)
        cla
        %Greg: updated inputs order
        [~] = PlotPhaseSpace([atc(ind1), iptg(ind1)], p,'LacImax',5e3,'TetRmax',3e3);
        for ind2 = 1:numcells
            plot(othercells(ind2).LacI(ind1),othercells(ind2).TetR(ind1),' *','markersize',10,'color','c')
        end
        plot(objective(1),objective(2),'w','Marker','pentagram')
        plot(LacI(ind1),TetR(ind1),' *b','markersize',10)
        drawnow
        %pause(.05)
    end
end

