function valvesplotscript(XP)

media = XP.media; % IPTG and ATC
% Plot the iptg and atc stairs as areas:
% [iptgstrX, iptgstrY] = stairs([0; media(:,3)]/3600,[media(:,1); media(end,1)]);
% [atcstrX, atcstrY] = stairs([0; media(:,3)]/3600,[media(:,2); media(end,2)]+125);
% h = area(iptgstrX, 100*iptgstrY, 0);
% set(h(1),'FaceColor',[0 1 1]);
% set(h(1),'EdgeColor',[0 1 1]);
% h = area(atcstrX, atcstrY, 125);
% set(h(1),'FaceColor',[1 0 1]);
% set(h(1),'EdgeColor',[1 0 1]);
hold on
stairs([0; media(:,3)]/3600,100*[media(:,1); media(end,1)],'color',[0 1 1],'linewidth',2);
stairs([0; media(:,3)]/3600,[media(:,2); media(end,2)]+125,'color',[1 0 1],'linewidth',2);

% Re-set position
set(gcf,'Position',[500 10 1000 1000]);
xlim([0 XP.timepoints(end)/3600])
ylim([0 225])
set(gca,'FontSize',12);

% Plot max and min lines:
line([xlim],[0     0],'Color','k')
line([xlim],[100 100],'Color','k')
line([xlim],[125 125],'Color','k')
line([xlim],[225 225],'Color','k')

% Cosmetics:
set(gca,'XTick',0:2:20);
set(gca,'YTick',[0 100 125 225]);
set(gca,'YTickLabel',{'0%' '1mM' '0' '100 ng/mL'});
set(gca,'Box','off');

xlabel('Time (h)')
ylabel('IPTG         aTC   ')