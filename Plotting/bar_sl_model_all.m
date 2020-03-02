% plot initMIP results

clear

%modscen='MIROC5_rcp85';
%modscen='initMIP';
modscen='M39_MIROC5_rcp85';

load ../Data/initMIP/ch_A5.mat
load(['../Data/initMIP/A5_iadsmb_trans_' modscen], 'prog');

ax=[0,100,-5,200];
%ylab='Sea-level contribution [mm]';
ylab='[mm]';

%colors=get(0,'DefaultAxesColorOrder');
colors = distinguishable_colors(36);

% color coded strings
igrpmodcol={};
for m=1:prog.n
    %\color[rgb]{0,0.5,0.5} text'
    igrpmodcol{m}=['\color[rgb]{',num2str(colors(m,1)),',',num2str(colors(m,2)),',',num2str(colors(m,3)),',0} ', prog.ids{ch.order(m)}];
end

% Plot init
figure(1)
set(gcf,'Position',[440   383   750   415])
%title('Sea-level contribution [SLE]');
hold on; box on
for m=1:prog.n;

    bar(m,prog.iadsmb(end,m),'FaceColor',colors(m,:));

end
hold off
% Make up
set(gca,'Xtick',1:1:prog.n)
set(gca,'XtickLabel',igrpmodcol,'XTickLabelRotation', 90)
set(gca,'FontSize', 16)
ylabel(ylab);
ax=axis;
axis([0 prog.n+1 ax(3) ax(4)])
ax=axis;
% save
print('-r300', '-dpng', ['Plots/A5_bar_trans_' modscen]);


% Plot remap
figure(2)
set(gcf,'Position',[440   383   750   415])
%title('Sea-level contribution [SLE]');
hold on; box on
for m=1:prog.n;

    bar(m,prog.iadsmb_re(end,m),'FaceColor',colors(m,:));

end
hold off
% Make up
set(gca,'Xtick',1:1:prog.n)
set(gca,'XtickLabel',igrpmodcol,'XTickLabelRotation', 90)
set(gca,'FontSize', 16)
ylabel(ylab);
axis(ax)
%ax=axis;
%axis([0 prog.n+1 ax(3) ax(4)])
% save
print('-r300', '-dpng', ['Plots/A5_bar_trans_remap_' modscen]);
posfull = get(gca, 'Position');

% Plot difference
figure(3)
set(gcf,'Position',[440   383   750   415])
%title('Sea-level contribution [SLE]');
hold on; box on
for m=1:prog.n;

    bar(m,prog.iadsmb(end,m)-prog.iadsmb_re(end,m),'FaceColor',colors(m,:));

end
hold off
% Make up
set(gca,'Xtick',1:1:prog.n)
set(gca,'XtickLabel',igrpmodcol,'XTickLabelRotation', 90)
set(gca,'FontSize', 16)
ylabel(ylab);
axis(ax)
axis([0 prog.n+1 ax(3) 70])
set(gca,'Position',[posfull(1),posfull(2),posfull(3),posfull(4)/(160/70)])
% save
print('-r300', '-dpng', ['Plots/A5_bar_trans_diff_' modscen]);

% Save results
dSLE = prog.iadsmb(end,:)-prog.iadsmb_re(end,:);
save(['../Data/initMIP/dSLE_' modscen '.mat'], 'dSLE');
