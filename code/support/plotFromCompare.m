function [] = plotFromCompare(predictions, data, name,fit)



yout = predictions.OutputData;
xout = data.OutputData;

ts = 1:length(yout);

figure; hold on; 
names = ["Actual"+name,...
        " prediction accuracy="+fit(1)];

    line1 = plot(ts, xout(:,1),'-r','linewidth',1.5);
    line2 = plot(ts, yout(:,1),':b','linewidth',1.5);
    set(gcf,'color','w');
    set(gca,'FontSize',14)
    ylabel("x_"+2)
    box on
% subplot(2,1,2); hold on;
%     line1 = plot(ts, xout(:,2),'-r','linewidth',1.5);
%     line2 = plot(ts, yout(:,2),':b','linewidth',1.5);
%     set(gcf,'color','w');
%     set(gca,'FontSize',14)
%     ylabel("x_"+3)
%     box on    
xlabel('Time Steps')
% Construct a Legend with the data from the sub-plots
hL = legend([line1,line2],names,'Orientation','horizontal');
% Programatically move the Legend
newPosition = [.52 .95 0 0];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);    
    
