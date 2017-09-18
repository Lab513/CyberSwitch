function DrawPresence(XP)
% This functions draws the presence heatmap of an experiment (i.e. parts of the state space that cells spend time in)

indxMstart = find(XP.timepoints>=XP.media(end-6,3),1,'first'); % Find beginning of last 2 periods
indxMstop = find(XP.timepoints>=XP.media(end-2,3),1,'first'); % Find end of last 2 periods
r = XP.rfp(indxMstart:indxMstop,:);
g = XP.gfp(indxMstart:indxMstop,:);
[n,c] = my2dhist([r(:),g(:)],100,100,150,100,[0 5e3],[0 2500]);
[percents] = percentify(n);
contourf(c(:,1), c(:,2), percents',[0 20 40 60 80 95 99.9999],'LineWidth',2); %99.9999 and not 100 because that would include the entire state space
colormap(flipud(parula))
xlabel('LacI-RFP')
ylabel('TetR-GFP')
hc = colorbar;
hc.Ticks = [ 20 40 60 80 95 99.9999];
hc.TickLabels = { '20' '40' '60' '80' '95' '100'};
ylabel(hc, 'Cells presence (%)')

function [percents] = percentify(n)

percents = zeros(size(n));
[descend, ~] = sort(n(:),1,'descend');

running_int = 0;

for ind1 = 1:numel(descend)
    mask = n <= descend(ind1);
    running_int = sum(sum(n(mask)))/sum(n(:));
    percents(mask) = running_int;
end

percents = 100*(1-percents);

function [N,C] = my2dhist(X,NBINS1,NBINS2,OVLP1,OVLP2,RNG1,RNG2)

C1 = linspace(RNG1(1),RNG1(2),NBINS1);
C2 = linspace(RNG2(1),RNG2(2),NBINS2);
X = X';

for ind1 = 1:numel(C1)
    cc1= C1(ind1);
    mi1 = cc1-OVLP1;
    ma1 = cc1+OVLP1;
    C(ind1,1) = cc1;
    for ind2 = 1:numel(C2)
        cc2= C2(ind2);
        mi2 = cc2-OVLP2;
        ma2 = cc2+OVLP2;
        xx = X(1,:) >= mi1 & X(1,:) < ma1 & X(2,:) >= mi2 & X(2,:) < ma2;
        N(ind1,ind2) = sum(double(xx(:)));
        C(ind2,2) = cc2;
    end
end