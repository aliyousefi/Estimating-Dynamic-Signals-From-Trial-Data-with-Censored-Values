function [Xk,Yn,Yb,In,Param] = ay_data_generator(N,OType,Range,a1,so,T)

%   Detailed explanation goes here
if nargin ==0
   %% Define data length
   N           = 100;  % Length of observartion
   a1          = 0.98; % AR(1) 
   Range       = 0.25;  % Range of Xh
   T           = 0.25; % Threshold
   OType       = 1;    % 1 means normal, 2 means Gamma
   so          = 0.001;
end

%% Set the Threshold
Param.T     = T;

%% Run the Hidden Model
Param.a1 = a1;
Param.a0 = Range *(1-Param.a1)*2;
Param.sh = Range^2 * (1-Param.a1*Param.a1);

% create the Xk
Xk    = zeros(N,1);
Xk(1) = Param.a0+sqrt(Param.sh) *randn;
for i=2:N
     Xk(i)=Param.a1*Xk(i-1)+Param.a0+sqrt(Param.sh) *randn;
end

%% Run the RT model
Yn = zeros(N,1);
% normal
if OType == 1
    Param.b1    = 1;
    Param.b0    = -0.6;
    Param.so    = so;
    
    % Generate data
    Yn = Param.b1 * Xk + Param.b0 + sqrt(Param.so) * randn(N,1);
end

% Gamma
if OType == 2
    Param.b1    = 1;
    Param.b0    = 0;
    Param.v     = 5;
    Param.alpha = 0;
   
    % Generate data
    Un  = exp(Param.b1 * Xk + Param.b0);
    Yn= Param.alpha+random('gamma',Param.v,Un/Param.v);
    In(find(Yn>T))=0;
end


%% Run the Binary model
Scale = 10.0;
%% Run the Binary model
BinaryExist = 1;
if BinaryExist
    Param.c2 = -1 * Scale;
    Param.c1 = 0.85 * Scale;
    Param.c0 = -0.25 * Scale;
    
    Yb  = zeros(N,1);
    
    Pb  = exp(Param.c2*Xk+Param.c1*exp(Yn)+Param.c0)./ (1+exp(Param.c2*Xk+Param.c1*exp(Yn)+Param.c0));
    Pt  = rand(N,1);
    ind = find(Pb>Pt);
    Yb(ind)= 1;
else
    Yb = [];
end

%% Find Threshold Index
Tn = length(T);
In = zeros(N,Tn);
for i=1:Tn
    ind = find(Yn<=T(i));
    In(ind,i) = 1;
end

end

