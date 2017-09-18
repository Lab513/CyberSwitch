%% This file can be used to re-plot most of the figures of the paper:

% First load up the parameters and add all the subpath of the directory:
% (You must be at the root of the directory)
addpath(genpath('.' ));
DATADIR = 'Data';
load('parameters2.mat');

% Then execute any of the following sections: (Again you must execute this script from the root of the repository.)
% I recommend executing this script section by section, some parts take time...

%% Figure 1

% C
% See main_param_search.m and plot_mean_optimised_results.m


% F
IndLvl = [20 .25]; % aTC (ng/ml), IPTG (mM)
figure('Name','Fig 1C','NumberTitle','off');
PlotPhaseSpace(IndLvl,p,'LacImax',3e3,'TetRmax',1.8e3, 'quiver', 'off', 'BKGD_image', 'none');
figure('Name','Fig 1C - insert','NumberTitle','off');
PlotPhaseSpace(IndLvl,p,'LacImin',200,'LacImax',1e3,'TetRmin',200,'TetRmax',500, 'quiver', 'on', 'BKGD_image', 'none');

%% Figure 2

% B
% To run a similar experiment, see mainControl.m
load(fullfile(DATADIR,'Simulated','PI_1','Data.mat'));
figure('Name','Fig 2B','NumberTitle','off');
subplot(4,1,1)
curvesplotscript(XP,'method','levels','which','controlled','LacImax',2100,'TetRmax',2100)
subplot(4,1,[2,3])
curvesplotscript(XP,'method','ratio','which','controlled')
subplot(4,1,4)
valvesplotscript(XP)


% C
load(fullfile(DATADIR,'Experimental','PI_1','Data.mat'));
figure('Name','Fig 2C','NumberTitle','off');
subplot(4,1,1)
curvesplotscript(XP,'method','levels','which','controlled','LacImax',2100,'TetRmax',2100)
subplot(4,1,[2,3])
curvesplotscript(XP,'method','ratio','which','controlled')
subplot(4,1,4)
valvesplotscript(XP)

% D
load(fullfile(DATADIR,'Experimental','PI_1','Data.mat'));
figure('Name','Fig 2D','NumberTitle','off');
curvesplotscript(XP,'method','color2D','which','controlled')

%% Figure 3

% B
load(fullfile(DATADIR,'Experimental','PI_1','Data.mat'));
figure('Name','Fig 3B','NumberTitle','off');
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)

% C
load(fullfile(DATADIR,'Experimental','BangBang_1','Data.mat'));
figure('Name','Fig 3C','NumberTitle','off');
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)

% E
% To run a similar experiment, see mainControl.m
load(fullfile(DATADIR,'Simulated','PI_1','Data.mat'));
figure('Name','Fig 3B','NumberTitle','off');
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)

% F
% To run a similar experiment, see simulateDynStim.m and the code for Fig
% 4B & 4G
load(fullfile(DATADIR,'Simulated','BangBang_1','Data.mat'));
figure('Name','Fig 3B','NumberTitle','off');
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)

%% Figure 4

% A
load(fullfile(DATADIR,'Experimental','DynStim_1','Data.mat'));
figure('Name','Fig 4A','NumberTitle','off');
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)

% B
% I didn't keep the simulation data for this one, but you can simulate the
% same experiment and get a similar result: (It can take some time)
XP = simulate_DynStim(p,[0 .5],[50 0],120,30,40,[20 .25],20*60,16);
figure('Name','Fig 4B','NumberTitle','off');
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)

% C
load(fullfile(DATADIR,'Experimental','DynStim_1','Data.mat'));
figure('Name','Fig 4C - 2D trajectory','NumberTitle','off');
curvesplotscript(XP,'method','color2D','which','1','LacImax',5000,'TetRmax',2500)
figure('Name','Fig 4C - presence overlay','NumberTitle','off');
DrawPresence(XP)


% D
IndLvl = [0 .5 120; 50 0 30]; % aTC (ng/ml), IPTG (mM), duration (min)
IndLvl_int = generate_integrative_inputs(p,IndLvl');
figure('Name','Fig 4D','NumberTitle','off');
PlotPhaseSpace(IndLvl_int',p,'LacImax',4e3,'TetRmax',1.5e3, 'quiver', 'off', 'BKGD_magnification_factor',40, 'BKGD_image', 'intensity', 'Nullclines', 'off');

% E
AllDynStim = dir(fullfile(DATADIR,'Experimental','DynStim_*'));
for ind1 = 1:numel(AllDynStim)
    load(fullfile(DATADIR,'Experimental',AllDynStim(ind1).name,'Data.mat'));
    indxMstart = find(XP.timepoints>=XP.media(end-6,3),1,'first'); % Find beginning of last 2 periods
    indxMstop = find(XP.timepoints>=XP.media(end-2,3),1,'first'); % Find end of last 2 periods
    rfpmean(ind1) = mean(mean(XP.rfp(indxMstart:indxMstop,:),2));
    gfpmean(ind1) = mean(mean(XP.gfp(indxMstart:indxMstop,:),2));
    ratio_comp = @(m) 0.01*sum(m(:,1).*m(:,3))/(sum(m(:,2).*m(:,3)));
    Durations = diff(XP.media(:,3));
    IPTGaTc(ind1) = ratio_comp([XP.media(2:3,1:2) Durations(1:2)]);
end
[IPTGaTc,indexes] = sort(IPTGaTc);
cm = getcolormap(log10(1./IPTGaTc),'cool');
figure('Name','Fig 4E','NumberTitle','off');
hold on;
for ind1 = 1:numel(IPTGaTc)
    plot(rfpmean(indexes(ind1)),gfpmean(indexes(ind1)),'Marker','.','Color',cm(ind1,:),'MarkerSize',20);
end
xlim([0 4e3]);
ylim([0 1.5e3]);
xlabel('LacI-RFP');
ylabel('TetR-GFP');

% F
figure('Name','Fig 4F','NumberTitle','off');
load(fullfile(DATADIR,'Experimental','DynStim_2','Data.mat'));
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)

% G
% I didn't keep the simulation data for this one, but you can simulate the
% same experiment and get a similar result:
XP = simulate_DynStim(p,[0 .5],[50 0],180,30,40,[20 .25],20*60,16);
figure('Name','Fig 4G','NumberTitle','off');
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)

%% Figure S3
% See main_param_search.m and plot_mean_optimised_results.m


%% Figure S4


%% Figure S5

% A
IndLvl = [0 1]; % aTC (ng/ml), IPTG (mM)
figure('Name','Supp Fig 5A','NumberTitle','off');
PlotPhaseSpace(IndLvl,p,'LacImax',4e3,'TetRmax',1.5e3, 'quiver', 'on', 'BKGD_image', 'none');

% B
IndLvl = [100 0]; % aTC (ng/ml), IPTG (mM)
figure('Name','Supp Fig 5B','NumberTitle','off');
PlotPhaseSpace(IndLvl,p,'LacImax',4e3,'TetRmax',1.5e3, 'quiver', 'on', 'BKGD_image', 'none');

% C
IndLvl = [20 .25]; % aTC (ng/ml), IPTG (mM)
figure('Name','Supp Fig 5C','NumberTitle','off');
PlotPhaseSpace(IndLvl,p,'LacImax',4e3,'TetRmax',1.5e3, 'quiver', 'on', 'BKGD_image', 'none');


%% Figure S6

% A, C-E (Use the )
[MSS, BSS, USS] = get_intersections_jb(p);
LSS = BSS(log(BSS(:,3)./BSS(:,4))> 1,:);
TSS = BSS(log(BSS(:,3)./BSS(:,4))<= 1,:);
figure('Name','Supp Fig 6 A,C,D,E','NumberTitle','off');
hold on
plot3(MSS(:,1),MSS(:,2),log10(MSS(:,3)./MSS(:,4)),' .b')
plot3(LSS(:,1),LSS(:,2),log10(LSS(:,3)./LSS(:,4)),' .r')
plot3(TSS(:,1),TSS(:,2),log10(TSS(:,3)./TSS(:,4)),' .g')
plot3(USS(:,1),USS(:,2),log10(USS(:,3)./USS(:,4)),' .', 'Color',[.7 .4 0])
xlabel('aTc')
ylabel('IPTG')
zlabel('LacI-RFP/TetR-GFP')

% B
[MSS_l, BSS_l, USS_l] = get_intersections_line(p);
LSS_l = BSS_l(log(BSS_l(:,3)./BSS_l(:,4))> 1,:);
TSS_l = BSS_l(log(BSS_l(:,3)./BSS_l(:,4))<= 1,:);
figure('Name','Supp Fig 6B','NumberTitle','off');
hold on
plot(MSS_l(:,1),log10(MSS_l(:,3)./MSS_l(:,4)),' .b')
plot(LSS_l(:,1),log10(LSS_l(:,3)./LSS_l(:,4)),' .r')
plot(TSS_l(:,1),log10(TSS_l(:,3)./TSS_l(:,4)),' .g')
plot(USS_l(:,1),log10(USS_l(:,3)./USS_l(:,4)),' .', 'Color',[.7 .4 0])
xlabel('aTc')
ylabel('LacI-RFP/TetR-GFP')


%% Figure S7

% B
load(fullfile(DATADIR,'Experimental','PI_2','Data.mat'));
figure('Name','Supp Fig 7B','NumberTitle','off');
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)

% C
load(fullfile(DATADIR,'Experimental','PI_2','Data.mat'));
figure('Name','Supp Fig 7C','NumberTitle','off');
curvesplotscript(XP,'method','color2D','which','controlled')

% E
load(fullfile(DATADIR,'Simulated','PI_2','Data.mat'));
figure('Name','Supp Fig 7E','NumberTitle','off');
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)

% F
load(fullfile(DATADIR,'Simulated','PI_2','Data.mat'));
figure('Name','Supp Fig 7F','NumberTitle','off');
curvesplotscript(XP,'method','color2D','which','controlled')

% H
load(fullfile(DATADIR,'Experimental','PI_3','Data.mat'));
figure('Name','Supp Fig 7H','NumberTitle','off');
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)

% I
load(fullfile(DATADIR,'Experimental','BangBang_2','Data.mat'));
figure('Name','Supp Fig 7I','NumberTitle','off');
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)

%% Figure S8

% A
load(fullfile(DATADIR,'Experimental','DynStim_3','Data.mat'));
figure('Name','Supp Fig 8A - Trajectory','NumberTitle','off');
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)
figure('Name','Supp Fig 8A - Presence','NumberTitle','off');
DrawPresence(XP)

% B
load(fullfile(DATADIR,'Experimental','DynStim_1','Data.mat'));
figure('Name','Supp Fig 8B - Trajectory','NumberTitle','off');
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)
figure('Name','Supp Fig 8B - Presence','NumberTitle','off');
DrawPresence(XP)

% C
load(fullfile(DATADIR,'Experimental','DynStim_2','Data.mat'));
figure('Name','Supp Fig 8C - Trajectory','NumberTitle','off');
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)
figure('Name','Supp Fig 8C - Presence','NumberTitle','off');
DrawPresence(XP)

% D
load(fullfile(DATADIR,'Experimental','DynStim_4','Data.mat'));
figure('Name','Supp Fig 8D - Trajectory','NumberTitle','off');
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)
figure('Name','Supp Fig 8D - Presence','NumberTitle','off');
DrawPresence(XP)

% E
load(fullfile(DATADIR,'Experimental','DynStim_5','Data.mat'));
figure('Name','Supp Fig 8E - Trajectory','NumberTitle','off');
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)
figure('Name','Supp Fig 8E - Presence','NumberTitle','off');
DrawPresence(XP)

%% Figure S9

load(fullfile(DATADIR,'Experimental','DynStimTooFast','Data.mat'));
figure('Name','Supp Fig 9','NumberTitle','off');
subplot(3,1,[1,2])
curvesplotscript(XP,'method','ratio','which','all')
subplot(3,1,3)
valvesplotscript(XP)