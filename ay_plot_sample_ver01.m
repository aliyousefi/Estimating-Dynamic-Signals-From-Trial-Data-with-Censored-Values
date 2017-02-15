close all
% generate data 
T= 0.9;
[Xk,Yn,Yb,In,Param] =ay_data_generator(100,1,0.25,0.95,0.08,log(T));
% sample range
Xs = -1:0.01:3;


figure(1)
plot(Xk,'b','LineWidth',3)
xlabel('Time')
ylabel('\it{X_k}');

%set(gca,'Color',[0.2 0.2 0.5]);

figure(2)
ind = find(In==0);
Yt = exp(Yn);
Yt(ind)=inf;
plot(Yt,'Color',[0 0 1],'LineWidth',3);
hold on;
plot(T*ones(100,1),'k--','LineWidth',0.3);
%Nums ={'','False', num2str(0.2),num2str(0.4),num2str(0.6),num2str(0.8),'True',''};
%set(gca,'YTick',linspace(-0.1,1.1,8),'YTickLabel',Nums);
%plot(ind,1,'kx','LineWidth',2)

ind = find(In==0);
Yb(ind)=inf;
plot(Yb,'o','Color',[0 0 1],'MarkerSize',6,'LineWidth',2);

ind = find(In==0);
plot(ind,T,'x','Color',[0 0.0 0],'LineWidth',2,'MarkerSize',8)
hold off;
xlabel('Time');
ylabel('Observation=\{{\it{Y_k} , \it{Z_k}}\}')
legend('Yk','T(sec)','Zk','Missed','Orientation','horizontal')

%set(gca,'Color',[0.2 0.2 0.5]);
ylim([-0.1 1.1])

save('keep_file','Xk','Yn','Yb');

% match to Xs
% ImagXk = zeros(100,1);
% for s=1:100
%     [temp,ind]= min(abs(Xk(s)-Xs));
%     ImagXk(s) = length(Xs)-ind;
% end

% Run the exact
figure(10)
Pxz = ay_smoothing(1,1,Yn,Yb,In(:,1),Xs,Param,log(T));

mse1 = ay_mse(Xk,Xs,Pxz);
[cnt1,lbnd,ubnd] = ay_count(Xk,Xs,Pxz);
        

MeanP = zeros(100,1);
for s=1:100
    mx = Xs * Pxz(s,:)';
    %[temp,ind]= min(abs(mx-Xs));
    MeanP(s) = mx; %length(Xs)-ind;
end
%MeanP = MeanP(end:-1:1);
%Pxz = Pxz(:,end:-1:1);
%h=imagesc(sqrt(Pxz'));
%colormap(flipud(bone))
%hold on

ay_plot_bound(1,1:100,MeanP,lbnd' ,ubnd' );hold on;
plot(Xk,'b','LineWidth',2);hold on;

num1 = min(lbnd-0.2);
num2 = max(ubnd+0.2);
ylim([num1 num2])

ind=find(In(:,1)==0);
plot(ind,num1+0.1,'x','MarkerSize',8,'LineWidth',2,'Color',[0 0.0 0]);

hold off;
% index = (linspace(Xs(num1),Xs(num2),6));
% Nums ={ num2str(index(1)), num2str(index(2)),num2str(index(3)),num2str(index(4)),num2str(index(5)),num2str(index(6))};
% set(gca,'YTick',linspace(num1,num2,6),'YTickLabel',Nums(end:-1:1));
xlabel('Time')
ylabel('\it{X_k}');
title('Exact Method','FontWeight','normal');
legend('95% HPD Rgn','\it{X_k_|_k}','\it{X_k}','Orientation','horizontal')




figure(11)
Pxz = ay_smoothing(1,2,Yn,Yb,In(:,1),Xs,Param,log(T));

mse2 = ay_mse(Xk,Xs,Pxz);
[cnt2,lbnd,ubnd] = ay_count(Xk,Xs,Pxz);
        

MeanP = zeros(100,1);
for s=1:100
    mx = Xs * Pxz(s,:)';
    %[temp,ind]= min(abs(mx-Xs));
    MeanP(s) = mx; %length(Xs)-ind;
end
%MeanP = MeanP(end:-1:1);
%Pxz = Pxz(:,end:-1:1);
%h=imagesc(sqrt(Pxz'));
%colormap(flipud(bone))
%hold on

ay_plot_bound(1,1:100,MeanP,lbnd' ,ubnd' );hold on;
plot(Xk,'b','LineWidth',2);hold on;

num1 = min(num1,min(lbnd-0.2));
num2 = max(num2,max(ubnd+0.2));


ind=find(In(:,1)==0);
plot(ind,num1+0.1,'x','MarkerSize',8,'LineWidth',2,'Color',[0 0.0 0]);

hold off;
% index = (linspace(Xs(num1),Xs(num2),6));
% Nums ={ num2str(index(1)), num2str(index(2)),num2str(index(3)),num2str(index(4)),num2str(index(5)),num2str(index(6))};
% set(gca,'YTick',linspace(num1,num2,6),'YTickLabel',Nums(end:-1:1));
xlabel('Time')
ylabel('\it{X_k}');
title('Data Imputation Method','FontWeight','normal');
legend('95% HPD Rgn','\it{X_k_|_k}','\it{X_k}','Orientation','horizontal')


figure(12)
Pxz = ay_smoothing(1,3,Yn,Yb,In(:,1),Xs,Param,log(T));

mse3 = ay_mse(Xk,Xs,Pxz);
[cnt3,lbnd,ubnd] = ay_count(Xk,Xs,Pxz);
        

MeanP = zeros(100,1);
for s=1:100
    mx = Xs * Pxz(s,:)';
    %[temp,ind]= min(abs(mx-Xs));
    MeanP(s) = mx; %length(Xs)-ind;
end
%MeanP = MeanP(end:-1:1);
%Pxz = Pxz(:,end:-1:1);
%h=imagesc(sqrt(Pxz'));
%colormap(flipud(bone))
%hold on

ay_plot_bound(1,1:100,MeanP,lbnd' ,ubnd' );hold on;
plot(Xk,'b','LineWidth',2);hold on;

num1 = min(num1,min(lbnd-0.2));
num2 = max(num2,max(ubnd+0.2));


ind=find(In(:,1)==0);
plot(ind,num1+0.1,'x','MarkerSize',8,'LineWidth',2,'Color',[0 0.0 0]);

hold off;
% index = (linspace(Xs(num1),Xs(num2),6));
% Nums ={ num2str(index(1)), num2str(index(2)),num2str(index(3)),num2str(index(4)),num2str(index(5)),num2str(index(6))};
% set(gca,'YTick',linspace(num1,num2,6),'YTickLabel',Nums(end:-1:1));
xlabel('Time')
ylabel('\it{X_k}');
title('Data Deletion Method','FontWeight','normal');
legend('95% HPD Rgn','\it{X_k_|_k}','\it{X_k}','Orientation','horizontal')

figure(13)
[temp,Mx,Sx] =ay_smoothing(1,4,Yn,Yb,In(:,1),Xs,Param,log(T));
for s=1:100
    Pxz(s,:)=exp(-(Xs-Mx(s)).^2/(2*Sx(s)));
    Pxz(s,:)=Pxz(s,:)/sum(Pxz(s,:));
end
mse4 = ay_mse(Xk,Xs,Pxz);
[cnt4,lbnd,ubnd] = ay_count(Xk,Xs,Pxz);
        

MeanP = zeros(100,1);
for s=1:100
    mx = Xs * Pxz(s,:)';
    %[temp,ind]= min(abs(mx-Xs));
    MeanP(s) = mx; %length(Xs)-ind;
end
%MeanP = MeanP(end:-1:1);
%Pxz = Pxz(:,end:-1:1);
%h=imagesc(sqrt(Pxz'));
%colormap(flipud(bone))
%hold on

ay_plot_bound(1,1:100,MeanP,lbnd' ,ubnd' );hold on;
plot(Xk,'b','LineWidth',2);hold on;

num1 = min(num1,min(lbnd-0.2));
num2 = max(num2,max(ubnd+0.2));


ind=find(In(:,1)==0);
plot(ind,num1+0.1,'x','MarkerSize',8,'LineWidth',2,'Color',[0 0.0 0]);

hold off;
% index = (linspace(Xs(num1),Xs(num2),6));
% Nums ={ num2str(index(1)), num2str(index(2)),num2str(index(3)),num2str(index(4)),num2str(index(5)),num2str(index(6))};
% set(gca,'YTick',linspace(num1,num2,6),'YTickLabel',Nums(end:-1:1));
xlabel('Time')
ylabel('\it{X_k}');
title('Gaussian Approximation Method','FontWeight','normal');
legend('95% HPD Rgn','\it{X_k_|_k}','\it{X_k}','Orientation','horizontal')

figure(1)
ylim([num1 num2])

figure(10)
ylim([num1 num2])

figure(11)
ylim([num1 num2])

figure(12)
ylim([num1 num2])

figure(13)
ylim([num1 num2])






% Run the exact
figure(20)
Pxz = ay_smoothing(1,1,Yn,[],In(:,1),Xs,Param,log(T));

mse1 = ay_mse(Xk,Xs,Pxz);
[cnt1,lbnd,ubnd] = ay_count(Xk,Xs,Pxz);
        

MeanP = zeros(100,1);
for s=1:100
    mx = Xs * Pxz(s,:)';
    %[temp,ind]= min(abs(mx-Xs));
    MeanP(s) = mx; %length(Xs)-ind;
end
%MeanP = MeanP(end:-1:1);
%Pxz = Pxz(:,end:-1:1);
%h=imagesc(sqrt(Pxz'));
%colormap(flipud(bone))
%hold on

ay_plot_bound(1,1:100,MeanP,lbnd' ,ubnd' );hold on;
plot(Xk,'b','LineWidth',2);hold on;

num1 = min(lbnd-0.2);
num2 = max(ubnd+0.2);
ylim([num1 num2])

ind=find(In(:,1)==0);
plot(ind,num1+0.1,'x','MarkerSize',8,'LineWidth',2,'Color',[0 0.0 0]);

hold off;
% index = (linspace(Xs(num1),Xs(num2),6));
% Nums ={ num2str(index(1)), num2str(index(2)),num2str(index(3)),num2str(index(4)),num2str(index(5)),num2str(index(6))};
% set(gca,'YTick',linspace(num1,num2,6),'YTickLabel',Nums(end:-1:1));
xlabel('Time')
ylabel('\it{X_k}');
title('Exact Method','FontWeight','normal');
legend('95% HPD Rgn','\it{X_k_|_k}','\it{X_k}','Orientation','horizontal')




figure(21)
Pxz = ay_smoothing(1,2,Yn,[],In(:,1),Xs,Param,log(T));

mse2 = ay_mse(Xk,Xs,Pxz);
[cnt2,lbnd,ubnd] = ay_count(Xk,Xs,Pxz);
        

MeanP = zeros(100,1);
for s=1:100
    mx = Xs * Pxz(s,:)';
    %[temp,ind]= min(abs(mx-Xs));
    MeanP(s) = mx; %length(Xs)-ind;
end
%MeanP = MeanP(end:-1:1);
%Pxz = Pxz(:,end:-1:1);
%h=imagesc(sqrt(Pxz'));
%colormap(flipud(bone))
%hold on

ay_plot_bound(1,1:100,MeanP,lbnd' ,ubnd' );hold on;
plot(Xk,'b','LineWidth',2);hold on;

num1 = min(num1,min(lbnd-0.2));
num2 = max(num2,max(ubnd+0.2));


ind=find(In(:,1)==0);
plot(ind,num1+0.1,'x','MarkerSize',8,'LineWidth',2,'Color',[0 0.0 0]);

hold off;
% index = (linspace(Xs(num1),Xs(num2),6));
% Nums ={ num2str(index(1)), num2str(index(2)),num2str(index(3)),num2str(index(4)),num2str(index(5)),num2str(index(6))};
% set(gca,'YTick',linspace(num1,num2,6),'YTickLabel',Nums(end:-1:1));
xlabel('Time')
ylabel('\it{X_k}');
title('Data Imputation Method','FontWeight','normal');
legend('95% HPD Rgn','\it{X_k_|_k}','\it{X_k}','Orientation','horizontal')


figure(22)
Pxz = ay_smoothing(1,3,Yn,[],In(:,1),Xs,Param,log(T));

mse3 = ay_mse(Xk,Xs,Pxz);
[cnt3,lbnd,ubnd] = ay_count(Xk,Xs,Pxz);
        

MeanP = zeros(100,1);
for s=1:100
    mx = Xs * Pxz(s,:)';
    %[temp,ind]= min(abs(mx-Xs));
    MeanP(s) = mx; %length(Xs)-ind;
end
%MeanP = MeanP(end:-1:1);
%Pxz = Pxz(:,end:-1:1);
%h=imagesc(sqrt(Pxz'));
%colormap(flipud(bone))
%hold on

ay_plot_bound(1,1:100,MeanP,lbnd' ,ubnd' );hold on;
plot(Xk,'b','LineWidth',2);hold on;

num1 = min(num1,min(lbnd-0.2));
num2 = max(num2,max(ubnd+0.2));


ind=find(In(:,1)==0);
plot(ind,num1+0.1,'x','MarkerSize',8,'LineWidth',2,'Color',[0 0.0 0]);

hold off;
% index = (linspace(Xs(num1),Xs(num2),6));
% Nums ={ num2str(index(1)), num2str(index(2)),num2str(index(3)),num2str(index(4)),num2str(index(5)),num2str(index(6))};
% set(gca,'YTick',linspace(num1,num2,6),'YTickLabel',Nums(end:-1:1));
xlabel('Time')
ylabel('\it{X_k}');
title('Data Deletion Method','FontWeight','normal');
legend('95% HPD Rgn','\it{X_k_|_k}','\it{X_k}','Orientation','horizontal')

figure(23)
[temp,Mx,Sx] =ay_smoothing(1,4,Yn,[],In(:,1),Xs,Param,log(T));
for s=1:100
    Pxz(s,:)=exp(-(Xs-Mx(s)).^2/(2*Sx(s)));
    Pxz(s,:)=Pxz(s,:)/sum(Pxz(s,:));
end
mse4 = ay_mse(Xk,Xs,Pxz);
[cnt4,lbnd,ubnd] = ay_count(Xk,Xs,Pxz);
        

MeanP = zeros(100,1);
for s=1:100
    mx = Xs * Pxz(s,:)';
    %[temp,ind]= min(abs(mx-Xs));
    MeanP(s) = mx; %length(Xs)-ind;
end
%MeanP = MeanP(end:-1:1);
%Pxz = Pxz(:,end:-1:1);
%h=imagesc(sqrt(Pxz'));
%colormap(flipud(bone))
%hold on

ay_plot_bound(1,1:100,MeanP,lbnd' ,ubnd' );hold on;
plot(Xk,'b','LineWidth',2);hold on;

num1 = min(num1,min(lbnd-0.2));
num2 = max(num2,max(ubnd+0.2));


ind=find(In(:,1)==0);
plot(ind,num1+0.1,'x','MarkerSize',8,'LineWidth',2,'Color',[0 0.0 0]);

hold off;
% index = (linspace(Xs(num1),Xs(num2),6));
% Nums ={ num2str(index(1)), num2str(index(2)),num2str(index(3)),num2str(index(4)),num2str(index(5)),num2str(index(6))};
% set(gca,'YTick',linspace(num1,num2,6),'YTickLabel',Nums(end:-1:1));
xlabel('Time')
ylabel('\it{X_k}');
title('Gaussian Approximation Method','FontWeight','normal');
legend('95% HPD Rgn','\it{X_k_|_k}','\it{X_k}','Orientation','horizontal')


figure(20)
ylim([num1 num2])

figure(21)
ylim([num1 num2])

figure(22)
ylim([num1 num2])

figure(23)
ylim([num1 num2])




stop=1;


% figure(2)
% ind = find(In==0);
% Yt  = exp(Yn);
% Yt(ind)=inf;
% Yb(ind)=inf;
% h=plotyy(1:100,Yt,1:100,Yb);
% set(h(1),'LineWidth',2);
% set(h(2),'LineStyleOrder','.');
% hold on;
% plot(ones(100,1),'r-','LineWidth',1)
% plot(ind,1,'kx','LineWidth',2)
% hold off;
% ind = find(In==1);
% plot(ind,1,'kx','LineWidth',2)


