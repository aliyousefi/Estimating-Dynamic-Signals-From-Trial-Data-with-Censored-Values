function [count,LowBound,UpBound,N] = ay_count(Xs,XPx,Px,Extra)

% check input arguments
if nargin <2
    disp('not enough nargin')
    return;
end

UpBound = zeros(length(Xs),1);
LowBound = zeros(length(Xs),1);
if nargin==3
    dx=XPx(2)-XPx(1);
    %% count
    count = 0;
    for i=1:length(Xs)
        Xn = ay_hpd(XPx,Px(i,:),0.95);
        xd = min(abs(Xn - Xs(i)));
        if xd <= dx 
            count = count + 1;
        end
        LowBound(i)=min(Xn);
        UpBound(i) =max(Xn);
    end
else  % 4 means normal
    count = 0;
    for i=1:length(Xs)
        Xmin = XPx(i)-2*sqrt(Px(i));
        Xmax = XPx(i)+2*sqrt(Px(i));
        if Xs(i)>= Xmin &&  Xs(i)<= Xmax
            count = count + 1;
        end
        LowBound(i)=Xmin;
        UpBound(i) =Xmax;
    end
    
end
N= length(Xs);