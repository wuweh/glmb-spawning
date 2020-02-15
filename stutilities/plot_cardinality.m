cardFig = figure;
t = 0:K;
maxCard = 10;%max(N_true);
stairs(1:K,N_true,'k-'); hold on; grid on; xlim([1,K]); ylim([0,maxCard+1]);
cardFig.Children.XTick = 10:10:K;
cardFig.Children.YTick = 0:2:maxCard+1;
set(gcf,'position',[1611         556         898         256]);

plot(1:K,hat_N,'k.','MarkerSize',8);
