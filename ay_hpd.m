function [XP] = ay_hpd(X,P,Prcnt)
% It retunrs points of HPD
% Input:
%     1. X, is distribution sample poitns
%     2. P, is the pdf - not necessarily to be normalized
%     3. Prcnt, is the HPD range
if nargin <2
    disp('not enough nargin')
    return;
end
if nargin ==2
    Prcnt = 0.95;
end

% Normalize P
P = P/(eps+sum(P));
% Sort P and then X
[Ps,ind]=sort(P,'descend');
Xs = X(ind);
% Find the Prcnt index on the sorted distribution
CPs = cumsum(Ps);
% Find CPs index close to HPD Precent
[temp,ind] = min(abs(Prcnt-CPs));
XP=Xs(1:ind);

