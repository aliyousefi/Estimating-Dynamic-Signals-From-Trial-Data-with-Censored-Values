T=linspace(log(0.4),log(2.4),18);

LEFT_BASE = 0.089226;

A  = load('normal_binary_count.txt');
Ns = A(:,6:6:end);
T  = 100-mean(Ns);

figure(1)
A  = load('normal.txt');
At =  sqrt(A(:,1:5:end));
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
%h=errorbar(T,mA,sA,'b','LineWidth',2,'MarkerSize',2);hold on;
h(1)=ay_plot_bound_a(1,T,mA,(mA-sA),(mA+sA));hold on;
plot(T,mA,'ko','LineWidth',1);



At =  sqrt(A(:,5:5:end));
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(2)=ay_plot_bound_a(2,T,mA,(mA-sA),(mA+sA));
plot(T,mA,'ko','LineWidth',1);


At =  sqrt(A(:,3:5:end));
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(3)=ay_plot_bound_a(3,T,mA,(mA-sA),(mA+sA));
plot(T,mA,'ko','LineWidth',1);


At =  sqrt(A(:,4:5:end));
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(4)=ay_plot_bound_a(4,T,mA,(mA-sA),(mA+sA));
plot(T,mA,'ko','LineWidth',1);hold off;


xlabel('Expected Percentage of Censored Data');
ylabel('RMSE');
legend([h(1),h(2),h(3),h(4)],'Exact','Imputation','Deletion','Approximate','Location','northwest');
title('Reaction Time Observation','FontWeight','Normal')
xlim([1 100])
grid

ylim([0.5 4]* LEFT_BASE)
h=gca;
set(h,'YTick',[0.5  1  1.5  2  2.5  3  3.5 4]* LEFT_BASE)

hold on
set(gca,'Box','off');  
axesPosition = get(gca,'Position'); 
hNewAxes = axes('Position',axesPosition,...  %# Place a new axes on top...
                'Color','none',...           %#   ... with no background color
                'YLim',[0.5 4],...            %#   ... and a different scale
                'YTick',[0.5  1  1.5  2  2.5  3  3.5 4],...
                'YAxisLocation','right',...  %#   ... located on the right
                'XTick',[],...               %#   ... with no x tick marks
                'Box','off');
ylabel(hNewAxes,'Normalized RMSE');
% Create second Y axes on the right.


figure(2)
A  = load('normal_count.txt');
At =  A(:,1:6:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(1)=ay_plot_bound_a(1,T,mA,(mA-sA),(mA+sA));hold on;
plot(T,mA,'ko','LineWidth',1);


At =  A(:,5:6:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(A(:,1))); 
h(2)=ay_plot_bound_a(2,T,mA,(mA-sA),(mA+sA));
plot(T,mA,'ko','LineWidth',1);


At =  A(:,3:6:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(3)=ay_plot_bound_a(3,T,mA,(mA-sA),(mA+sA));
plot(T,mA,'ko','LineWidth',1);


At =  A(:,4:6:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(4)=ay_plot_bound_a(4,T,mA,(mA-sA),(mA+sA));
plot(T,mA,'ko','LineWidth',1);


xlabel('Expected Percentage of Censored Data');
ylabel('Coverage of 95% HPD Region');
legend([h(1),h(2),h(3),h(4)],'Exact','Imputation','Deletion','Approximate','Location','southwest');
%title('Reaction Time Observation','FontWeight','Normal')
%title('Reaction Time Observation')
xlim([1 100])
ylim([40 100])
grid on


% figure(13);
% At =  A(:,5:5:end);
% mA= mean(At);
% sA= std(At); sA = sA/sqrt(length(At(:,1))); 
% errorbar(T,mA,sA,'b','LineWidth',2);
% xlabel('Number of observed data');
% ylabel('Number of Observed Samples');
% 
% 
% 



figure(3)
A  = load('normal_binary.txt');
At =  sqrt(A(:,1:5:end));
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
%h=errorbar(T,mA,sA,'b','LineWidth',2,'MarkerSize',2);hold on;
h(1)=ay_plot_bound_a(1,T,mA,(mA-sA),(mA+sA));hold on;
plot(T,mA,'ko','LineWidth',1);



At =  sqrt(A(:,5:5:end));
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(2)=ay_plot_bound_a(2,T,mA,(mA-sA),(mA+sA));
plot(T,mA,'ko','LineWidth',1);


At =  sqrt(A(:,3:5:end));
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(3)=ay_plot_bound_a(3,T,mA,(mA-sA),(mA+sA));
plot(T,mA,'ko','LineWidth',1);


At =  sqrt(A(:,4:5:end));
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(4)=ay_plot_bound_a(4,T,mA,(mA-sA),(mA+sA));
plot(T,mA,'ko','LineWidth',1);hold off;

xlabel('Expected Percentage of Censored Data');
ylabel('RMSE');
legend([h(1),h(2),h(3),h(4)],'Exact','Imputation','Deletion','Approximate','Location','northwest');
title('Reaction Time plus Binary Decision Observation','FontWeight','Normal')
xlim([1 100])
grid

ylim([0.5 4]* LEFT_BASE)
h=gca;
set(h,'YTick',[0.5  1  1.5  2  2.5  3  3.5 4]* LEFT_BASE)

hold on
set(gca,'Box','off');  
axesPosition = get(gca,'Position'); 
hNewAxes = axes('Position',axesPosition,...  %# Place a new axes on top...
                'Color','none',...           %#   ... with no background color
                'YLim',[0.5 4],...            %#   ... and a different scale
                'YTick',[0.5  1  1.5  2  2.5  3  3.5 4],...
                'YAxisLocation','right',...  %#   ... located on the right
                'XTick',[],...               %#   ... with no x tick marks
                'Box','off');
ylabel(hNewAxes,'Normalized RMSE');
% Create second Y axes on the right.


figure(4)
A  = load('normal_binary_count.txt');
At =  A(:,1:6:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(1)=ay_plot_bound_a(1,T,mA,(mA-sA),(mA+sA));hold on;
plot(T,mA,'ko','LineWidth',1);


At =  A(:,5:6:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(A(:,1))); 
h(2)=ay_plot_bound_a(2,T,mA,(mA-sA),(mA+sA));
plot(T,mA,'ko','LineWidth',1);


At =  A(:,3:6:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(3)=ay_plot_bound_a(3,T,mA,(mA-sA),(mA+sA));
plot(T,mA,'ko','LineWidth',1);


At =  A(:,4:6:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(4)=ay_plot_bound_a(4,T,mA,(mA-sA),(mA+sA));
plot(T,mA,'ko','LineWidth',1);


xlabel('Expected Percentage of Censored Data');
ylabel('Coverage of 95% HPD Region');
legend([h(1),h(2),h(3),h(4)],'Exact','Imputation','Deletion','Approximate','Location','southwest');
%title('Reaction Time plus Binary Decision Observation','FontWeight','Normal')
%title('Reaction Time Observation')
xlim([1 100])
ylim([40 100])
grid on
