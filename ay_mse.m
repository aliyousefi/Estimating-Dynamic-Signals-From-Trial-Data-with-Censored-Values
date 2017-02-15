function [MSE] = ay_mse(Xc,Xs,Px)
% This function calculates the Mean Squared Error between Prediction and
% Correct one
% Input:
%     1.Xc, is the unobserved X
%     1.Xs, is the distribution sample points  

% check input arguments
if nargin <2
    disp('not enough nargin')
    return;
end

%% run the MSE
mse = zeros(length(Xc),1);
for i=1:length(Xc)
%     Xt = Xs - Xc(i);
%     Xt = Xt.*Xt;
%     mse(i)=sum(Xt.*Px(i,:));
    
    mx = sum(Xs.*Px(i,:));   
    mse(i)= (Xc(i)-mx)^2;

end
MSE = mean(mse);