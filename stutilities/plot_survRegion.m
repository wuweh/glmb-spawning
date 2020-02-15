XYfig = figure; 
% set(gcf,'position',[2511           1         690         661]);
set(gcf,'position',[2639         296         566         518]);

orng = [255, 128, 0]./255;

%--- range of the plots
limit2= 1000*[ -1 1 -1 1 ];
ms1 = 9;

LS = {'k-','k-','k-'};
LW = 1.2;
for tidx = 1:scenario.trackCount % loop over tracks
    itrack = scenario.truthPlot.X{tidx};
    pos_itrack = model.C_posn*reshape(itrack(~isnan(itrack)),x_dim,[]);
    if ismember(tidx,[1,4,10])
        ls = LS{1};
        hleg(1) =  plot(pos_itrack(1,:),pos_itrack(2,:),ls,'Linewidth',LW); hold on; grid on; axis(limit2); axis equal; xlim([-1000,1000]); ylim([-1000,1000]);
    
    elseif ismember(tidx,[2,5,11])
        ls = LS{2};
        hleg(2) =  plot(pos_itrack(1,:),pos_itrack(2,:),ls,'Linewidth',LW); hold on; grid on; axis(limit2); axis equal; xlim([-1000,1000]); ylim([-1000,1000]);
    
    elseif ismember(tidx,[3,6,12])
        ls = LS{3};
        hleg(3) =  plot(pos_itrack(1,:),pos_itrack(2,:),ls,'Linewidth',LW); hold on; grid on; axis(limit2); axis equal; xlim([-1000,1000]); ylim([-1000,1000]);
    
    else
        ls = 'k-';
        hleg(4) =  plot(pos_itrack(1,:),pos_itrack(2,:),ls,'Linewidth',LW); hold on; grid on; axis(limit2); axis equal; xlim([-1000,1000]); ylim([-1000,1000]);
    end
    plot(pos_itrack(1,:),pos_itrack(2,:),ls,'Linewidth',LW); hold on; grid on; axis(limit2); axis equal; xlim([-1000,1000]); ylim([-1000,1000]);
    if isequal(scenario.truthPlot.type{tidx},'birth')
        plot(pos_itrack(1,1),pos_itrack(2,1),'ko','MarkerSize',ms1,'Linewidth',LW);
    else
        plot(pos_itrack(1,1),pos_itrack(2,1),'ks','MarkerSize',ms1,'Linewidth',LW);
    end
    plot(pos_itrack(1,end),pos_itrack(2,end),'k^','MarkerSize',ms1,'Linewidth',LW);
end % loop over tracks
% legend(hleg,'Birth Region 1','Birth Region 2','Birth Region 3','No Offspring','Location','SouthEast');

for idx = 1:length(connect)
    xc = connect{idx}(:,1);
    yc = connect{idx}(:,2);
    plot(xc,yc,'k:','LineWidth',1.6);
    
end

h(1) = plot(-2000,-2000,'k-','LineWidth',LW);
h(2) = plot(-2000,-2000,'k:','LineWidth',1.6);

H = legend(h,'Target Track','Parent - Spawn Lineage');
set(H,'FontSize',12);


grid on
% set(gcf,'paperposition',[0,0,8,8]);
saveas(gcf,'scenario','fig');

% Pspawn = model.C_posn*Xspawn;
% plot(Pspawn(1,:),Pspawn(2,:),'k-');
% Pspawn = model.C_posn*Xs_tmp
% plot(Pspawn(1,:),Pspawn(2,:),'k-');


xlabel('x coordinate (m)','fontsize',15); ylabel('y coordinate (m)','fontsize',15);
axis(limit2); axis square;
% axis([-462.546816479401                       500         -473.708938761691          495.974854112924]);



% %% Plot true cardinality
% cardFig = figure;
% t = 0:K;
% maxCard = 10;%max(N_true);
% stairs(1:K,N_true,'k-'); hold on; grid on; xlim([1,K]); ylim([0,maxCard+1]);
% cardFig.Children.XTick = 10:10:K;
% cardFig.Children.YTick = 0:2:maxCard+1;
% set(gcf,'position',[1611         556         898         256]);

% %% Set up OSPA plot
% ospaFig = figure;
% t = 0:K;
% stairs(-100,-100,'w.'); hold on; grid on; xlim([1,K]); ylim([0,100]);
% ospaFig.Children.XTick = 10:10:K;
% ospaFig.Children.YTick = 0:10:100;
% set(gcf,'position',[1612         214         898         256]);
