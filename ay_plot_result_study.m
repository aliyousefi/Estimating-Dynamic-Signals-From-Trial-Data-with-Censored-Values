T=linspace(log(0.2),log(2),18);

A  = load('normal_binary_count.txt');
Ns = A(:,5:5:end);
T  = 100-mean(Ns);

figure(1)
A  = load('normal.txt');
At =  A(:,1:4:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
%h=errorbar(T,mA,sA,'b','LineWidth',2,'MarkerSize',2);hold on;
h(1)=ay_plot_bound_a(1,T,mA,(mA-sA),(mA+sA));hold on;




At =  A(:,2:4:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(2)=ay_plot_bound_a(2,T,mA,(mA-sA),(mA+sA));



At =  A(:,3:4:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(3)=ay_plot_bound_a(3,T,mA,(mA-sA),(mA+sA));



At =  A(:,4:4:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(4)=ay_plot_bound_a(4,T,mA,(mA-sA),(mA+sA));hold off;

xlabel('Expected Percentage of Missing Data');
ylabel('RMSE');
legend([h(1),h(2),h(3),h(4)],'Exact','Imputation','Deletion','Approximate');
title('Reaction Time Observation','FontWeight','Normal')
xlim([1 100])
grid


figure(2)
A  = load('normal_count.txt');
At =  A(:,1:5:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(1)=ay_plot_bound_a(1,T,mA,(mA-sA),(mA+sA));hold on;


At =  A(:,2:5:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(A(:,1))); 
h(2)=ay_plot_bound_a(2,T,mA,(mA-sA),(mA+sA));


At =  A(:,3:5:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(3)=ay_plot_bound_a(3,T,mA,(mA-sA),(mA+sA));


At =  A(:,4:5:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(4)=ay_plot_bound_a(4,T,mA,(mA-sA),(mA+sA));


xlabel('Expected Percentage of Missing Data');
ylabel('Coverage of 95% HPD Region');
legend([h(1),h(2),h(3),h(4)],'Exact','Imputation','Deletion','Approximate');
%title('Reaction Time Observation','FontWeight','Normal')
%title('Reaction Time Observation')
xlim([1 100])
grid


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
At =  A(:,1:4:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
%h=errorbar(T,mA,sA,'b','LineWidth',2,'MarkerSize',2);hold on;
h(1)=ay_plot_bound_a(1,T,mA,(mA-sA),(mA+sA));hold on;




At =  A(:,2:4:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(2)=ay_plot_bound_a(2,T,mA,(mA-sA),(mA+sA));



At =  A(:,3:4:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(3)=ay_plot_bound_a(3,T,mA,(mA-sA),(mA+sA));



At =  A(:,4:4:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(4)=ay_plot_bound_a(4,T,mA,(mA-sA),(mA+sA));hold off;

xlabel('Expected Percentage of Missing Data');
ylabel('RMSE');
legend([h(1),h(2),h(3),h(4)],'Exact','Imputation','Deletion','Approximate');
title('Reaction Time plus Binary Decision Observation','FontWeight','Normal')
xlim([1 100])
grid


figure(4)
A  = load('normal_binary_count.txt');
At =  A(:,1:5:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(1)=ay_plot_bound_a(1,T,mA,(mA-sA),(mA+sA));hold on;


At =  A(:,2:5:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(A(:,1))); 
h(2)=ay_plot_bound_a(2,T,mA,(mA-sA),(mA+sA));


At =  A(:,3:5:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(3)=ay_plot_bound_a(3,T,mA,(mA-sA),(mA+sA));


At =  A(:,4:5:end);
mA= mean(At);
sA= std(At); sA = 2*sA/sqrt(length(At(:,1))); 
h(4)=ay_plot_bound_a(4,T,mA,(mA-sA),(mA+sA));


xlabel('Expected Percentage of Missing Data');
ylabel('Coverage of 95% HPD Region');
legend([h(1),h(2),h(3),h(4)],'Exact','Imputation','Deletion','Approximate');
%title('Reaction Time plus Binary Decision Observation','FontWeight','Normal')
%title('Reaction Time Observation')
xlim([1 100])
grid
