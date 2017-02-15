function ay_simulator()

%% define data length
N           = 100;  % Length of observartion
BinaryExist = 1;    % 1 means binary exist
Range       = 0.5;  % Range of Xh
T           = 0.25;  % Threshold
OType       = 1;    % 1 means normal, 2 means Gamma
Xs          = -5*Range:0.1:5*Range;

%% Set the Threshold
Param.T     = T;


%% Run the Hidden Model
Param.a1 = 0.98;
Param.a0 = 0;
Param.sh = Range^2 * (1-Param.a1*Param.a1);

% create the Xk
Xk    = zeros(N,1);
Xk(1) = sqrt(Param.sh) *randn;
for i=2:N
     Xk(i)=Param.a1*Xk(i-1)+Param.a0+sqrt(Param.sh) *randn;
end

%% Run the RT model
Yn = zeros(N,1);
In = ones(N,1);
% normal
if OType == 1
    Param.b1    = 0.8;
    Param.b0    = 0;
    Param.so    = 0.001;
    
    % Generate data
    Yn = Param.b1 * Xk + Param.b0 + sqrt(Param.so) * randn(N,1);
    In(find(Yn>T))=0;
end
% Gamma
if OType == 2
    Param.b1    = 0.8;
    Param.b0    = 0;
    Param.v     = 5;
    Param.alpha = 0;
   
    % Generate data
    Un  = exp(Param.b1 * Xk + Param.b0);
    Yn= Param.alpha+random('gamma',Param.v,Un/Param.v);
    In(find(Yn>T))=0;
end


%% Run the Binary model
if BinaryExist
    Param.c1 = -1;
    Param.c0 = 0;
    
    Yb  = zeros(N,1);
    
    Pb  = exp(Param.c1*Xk+Param.c0)./ (1+exp(Param.c1*Xk+Param.c0));
    Pt  = rand(N,1);
    ind = find(Pb>Pt);
    Yb(ind)= 1;
else
    Yb = [];
end


%% Run the algorithm
CType = 5;
[temp1,temp2,Mix_FMx,Mix_FLx,Mix_FUx,Mix_SMx,Mix_SLx,Mix_SUx] = ay_smoothing(OType,CType,Yn,Yb,In,Xs,Param);

CType = 4;
[temp1,temp2,Gauss_FMx,Gauss_FLx,Gauss_FUx,Gauss_SMx,Gauss_SLx,Gauss_SUx] = ay_smoothing(OType,CType,Yn,Yb,In,Xs,Param);

CType = 1;
[Exact_Qxz,Exact_Pxz,Exact_FMx,Exact_FLx,Exact_FUx,Exact_SMx,Exact_SLx,Exact_SUx] = ay_smoothing(OType,CType,Yn,Yb,In,Xs,Param);

CType = 2;
[Impute_Qxz,Impute_Pxz,Impute_FMx,Impute_FLx,Impute_FUx,Impute_SMx,Impute_SLx,Impute_SUx,Impute_Yn,Impute_Yb] = ay_smoothing(OType,CType,Yn,Yb,In,Xs,Param);

CType = 3;
[Drop_Qxz,Drop_Pxz,Drop_FMx,Drop_FLx,Drop_FUx,Drop_SMx,Drop_SLx,Drop_SUx] = ay_smoothing(OType,CType,Yn,Yb,In,Xs,Param);




% Figure
figure(1)
subplot(2,1,1)
plot(Yn);
hold on;
plot(Yb,'o');
ind=find(In==0);
plot(ind,In(ind)-0.2,'*');
legend('Yn','Biary','Missing')
plot(1:N,T*ones(N,1));
hold off;

subplot(2,1,2)
plot(Xk,'LineWidth',2);hold on;
plot(Exact_SMx,'LineWidth',2);
plot(Impute_SMx,'LineWidth',2);
plot(Drop_SMx,'LineWidth',2);
plot(Gauss_SMx,'LineWidth',2);
plot(Mix_SMx,'LineWidth',2);
legend('Xk','Exact','Impute','Drop','Gauss','Mix');
ind=find(In==0);
plot(ind,In(ind),'*');
hold off;
title('Xk')

figure(2)
plot(Yn);hold on;
plot(Impute_Yn);
legend('Yn','Imput-Yn');
ind=find(In==0);
plot(ind,In(ind),'*');
hold off;
title('Impute ')
hold off;







