function [Pxz,FMx,FSx,Yn,Yb] = ay_smoothing(OType,CType,Yn,Yb,In,Xs,Param,Ts)
%% Exact Filtering
% OType : observation type "Normal" or "Gamma"
% CType : censor type "Exact","Impute","Drop", or "Gaussian Approximate"
% Yn is the continuous observation
% Yb is the bernoulli observation
% In is the indicator function for the missing data
% Xs is the data samples
% Param is the model parameters
IMPUT_MAX = 10;
Pxz = [];
FMx = [];
FSx = [];

if nargin==0
    % Obseravtiont type
    OType = 1;  % normal (1), Gamma (2)
    % Data censor type
    CType = 2;  % exact (1), impute (2), and drop (3), and Gaussian approximate (4) 
        
    % Parameters for the hidden state space model
    Param.a1    = 0.98;
    Param.a0    = 0.02;
    Param.sh    = 0.05;
    % Paameters of the observation
    if OType == 1
        Param.b1 = 0.5;
        Param.b0 = 0;
        Param.so = 0.1;
    else
        Param.b1 = 0.5;
        Param.b0 = 0;
        Param.v  = 1;
        Param.alpha  = 0.1;
    end
    % Bernoulli part
    Param.c1    = 1;
    Param.c2    = -1;
    Param.c0    =  0;
    % Threshold level
    Param.T     = 0.5;
    % Continuous input
    Yn    = rand(100,1);
    % Binary input
    Yb    = zeros(100,1);
    % Censored data index
    In    = zeros(100,1);
    % Data Samples - possible input range
    Xs    = -3:0.1:3;
end

%% Assign model paramaters
% Parameters of the hidden stae
a1    = Param.a1;
a0    = Param.a0;
sh    = Param.sh;
% Paameters of the observation
if OType==1
    b1    = Param.b1;
    b0    = Param.b0;
    so    = Param.so;
else
    b1    = Param.b1;
    b0    = Param.b0;
    v     = Param.v;
    alpha = Param.alpha;
end
% Bernoulli part
ExistYb = 0;
if isempty(Yb)==0
    ExistYb = 1;
    c1    = Param.c1;
    c2    = Param.c2;
    c0    = Param.c0;
end
% Threshold level
T     = Ts;

if nargin==0
    % create the Xk
    Xk = zeros(100,1);
    Xk(1) = sqrt(sh) *randn;
    for i=2:100
        Xk(i)=a1*Xk(i-1)+a0+sqrt(sh) *randn;
    end
    
    if OType == 1
        % Yn
        Yn = b1 * Xk + b0 + sqrt(so) *randn(100,1);
    else
        Un = exp(b1*Xk+b0);
        Yn=alpha+random('gamma',v,Un/v);
    end
    
    % Yb
    if isempty(Yb)==0
        if OType == 1
            Pb  = exp(c2*Xk+c1*exp(Yn)+c0)./ (1+exp(c2*Xk+c1*exp(Yn)+c0));
        else
            Pb  = exp(c2*Xk+c1*exp(Yn)+c0)./ (1+exp(c2*Xk+c1*exp(Yn)+c0));
        end
        Pt  = rand(100,1);
        ind = find(Pb>Pt);
        Yb(ind)= 1;
    end
    
    % In
    In  = ones(100,1);
    ind = find(Yn>T);
    In(ind)= 0;
    
    %% This is for simulation
    YnBackup = Yn;

end

N = length(Yn);

if CType~=4    % exact, imput and ignore methods
    %% Filtering section: calculate the exact posterior
    Pxz= zeros(N,length(Xs));
    for i=1:N
        if In(i)==0 && CType == 2   % Impute
            if OType == 1
                % calculate one step prediction
                Px = zeros(length(Xs),1);
                for t=1:length(Xs)
                    if i == 1
                        Px(t) = pdf('normal',Xs(t),a0,sqrt(sh));
                    else
                        Pa    = pdf('normal',Xs(t),a1*Xs+a0,sqrt(sh));
                        Px(t) = Pa* Pxz(i-1,:)';
                    end
                end
                CPx = cumsum(Px);
                Px = Px / sum(Px);
                % Sample from the Px, find the first sample for
                % which the result is above T
                [temp,ui] = min(abs(rand-CPx));
                ur = Xs(ui);
                
                ms  =  b1*ur + b0;
                vs  =  sqrt(so);
                ds  = T:0.0001:(T+10*vs);
                gpx = cumsum(normpdf(ds,ms,vs))/sum(normpdf(ds,ms,vs));
                [temp,ti] = min(abs(rand-gpx));
                or = ds(ti);
                % create the sample
                Yn(i)= or;
                if ExistYb
                    pt   = exp(c2*ur+c1*exp(Yn(i))+c0)./(1+exp(c2*ur+c1*exp(Yn(i))+c0));
                    Yb(i)= 0;
                    if pt < rand
                        Yb(i) = 1;
                    end
                end
            end
        end
        if In(i)==0 && CType == 5   % Multiple Impute
                if OType == 1
                    % calculate one step prediction
                    Px = zeros(length(Xs),1);
                    for t=1:length(Xs)
                        if i == 1
                            Px(t) = pdf('normal',Xs(t),a0,sqrt(sh));
                        else
                            Pa    = pdf('normal',Xs(t),a1*Xs+a0,sqrt(sh));
                            Px(t) = Pa* Pxz(i-1,:)';
                        end
                    end
                    CPx = cumsum(Px);
                    Px = Px / sum(Px);
                    % Sample from the Px, find the first sample for
                    % which the result is above T
                    samp = 1;
                    while samp <= IMPUT_MAX
                       [temp,ui] = min(abs(rand-CPx));
                       ur = Xs(ui);
                       % ur is the sample of X
                       % now, generate Y sample
                       % ur is the sample of X
                       % now, generate Y sample
                       ms  =  b1*ur + b0;
                       vs  =  sqrt(so);
                       ds  = T:0.0001:(T+10*vs);
                       gpx = cumsum(normpdf(ds,ms,vs))/sum(normpdf(ds,ms,vs));
                       [temp,ti] = min(abs(rand-gpx));
                       or = ds(ti);
                       tYn(samp)= or;
                       % create the sample
                       if ExistYb
                            pt   = exp(c2*ur+c1*exp(Yn(i))+c0)./(1+exp(c2*ur+c1*exp(Yn(i))+c0));
                            tYb(samp)= 0;
                            if pt < rand
                                tYb(samp) = 1;
                            end
                       end
                       samp = samp+1;
                    end
                end
            if OType == 2
                % calculate one step prediction
                Px = zeros(length(Xs),1);
                for t=1:length(Xs)
                    if i == 1
                        Px(t) = pdf('normal',Xs(t),a0,sqrt(sh));
                    else
                        Pa    = pdf('normal',Xs(t),a1*Xs+a0,sqrt(sh));
                        Px(t) = Pa* Pxz(i-1,:)';
                    end
                end
                CPx = cumsum(Px);
                Px = Px / sum(Px);
                
                % Sample from the Px, find the first sample for
                % which the result is above T
                samp = 1;
                while samp
                        [temp,ui] = min(abs(rand-CPx));
                        ur = Xs(ui);
                        % ur is the sample of X
                        % now, generate Y sample
                        ut = exp(b1*ur+b0);
                        or = random('gamma',v,ut/v);
                        if or+alpha > T
                            samp = 0;
                        end
                end
                Yn(i)= or+alpha;
                % run the filtering
                if ExistYb
                    pt = exp(c2*ur+c1*exp(Yn(i))+c0)./(1+exp(c2*ur+c1*exp(Yn(i))+c0));
                    Yb(i)= 0;
                    if pt > rand
                        Yb(i)=1;
                    end
                end
            end
        end
        %% Main body
        if OType == 1   % Normal
            if In(i)==1 % Observed
                if ExistYb
                    pt = exp(c2*Xs+c1*exp(Yn(i))+c0)./(1+exp(c2*Xs+c1*exp(Yn(i))+c0));
                    p1 = pdf('normal',Yn(i),b1*Xs+b0,sqrt(so)).*(pt.^Yb(i)).*((1-pt).^(1-Yb(i)));
                else
                    p1 = pdf('normal',Yn(i),b1*Xs+b0,sqrt(so));
                end
            else        % Censored data
                if CType == 1 % exact
                    p1 = cdf('normal',T,b1*Xs+b0,sqrt(so),'upper');
                end
                if CType == 3 % ignore
                    p1 = ones(size(Xs));
                end    
                if CType == 2 % sampling (imputation)
                    if ExistYb
                        % run the filtering
                        pt = exp(c2*Xs+c1*exp(Yn(i))+c0)./(1+exp(c2*Xs+c1*exp(Yn(i))+c0));
                        p1 = pdf('normal',Yn(i),b1*Xs+b0,sqrt(so)).*(pt.^Yb(i)).*((1-pt).^(1-Yb(i)));
                    else
                        p1 = pdf('normal',Yn(i),b1*Xs+b0,sqrt(so));
                    end
                end
                if CType == 5 % sampling (imputation)
                    if ExistYb
                        % run the filtering
                        pt = exp(c2*Xs+c1*exp(tYn(1))+c0)./(1+exp(c2*Xs+c1*exp(tYn(1))+c0));
                        p1 = pdf('normal',tYn(1),b1*Xs+b0,sqrt(so)).*(pt.^tYb(1)).*((1-pt).^(1-tYb(1)));
                    else
                        p1 = pdf('normal',tYn(1),b1*Xs+b0,sqrt(so));
                    end
                    for ll=2:IMPUT_MAX
                        if ExistYb
                            % run the filtering
                            pt = exp(c2*Xs+c1*exp(tYn(ll))+c0)./(1+exp(c2*Xs+c1*exp(tYn(ll))+c0));
                            p1 = p1+pdf('normal',tYn(ll),b1*Xs+b0,sqrt(so)).*(pt.^tYb(ll)).*((1-pt).^(1-tYb(ll)));
                        else
                            p1 = p1+pdf('normal',tYn(ll),b1*Xs+b0,sqrt(so));
                        end
                    end
                    p1 = p1/IMPUT_MAX;
                end
            end
        end
        if OType == 2   % Gamma
            if In(i)==1 % Observed
                    if ExistYb
                        pt = exp(c2*Xs+c1*exp(Yn(i))+c0)./(1+exp(c2*Xs+c1*exp(Yn(i))+c0));
                        u  = exp(b1*Xs+b0);
                        p1 = gamma_pdf(Yn(i)-alpha,u,v).*(pt.^Yb(i)).*((1-pt).^(1-Yb(i)));
                    else
                        p1 = gamma_pdf(Yn(i)-alpha,u,v);
                    end
            else        % Censored 
                    if CType == 1 % exact
                        u  = exp(b1*Xs+b0);
                        p1 = gammainc(v*(T-alpha)/u,v,'upper');
                    end
                    if CType == 2 % sampling
                        % run the filtering
                        if ExistYb
                            pt = exp(c2*Xs+c1*exp(Yn(i))+c0)./(1+exp(c2*Xs+c1*exp(Yn(i))+c0));
                            u  = exp(b1*Xs+b0);
                            p1 = gamma_pdf(Yn(i)-alpha,u,v).*(pt^Yb(i)).*((1-pt).^(1-Yb(i)));
                        else
                            p1 = gamma_pdf(Yn(i)-alpha,u,v);
                        end
                    end
                    if CType == 3 % ignor
                        p1 = ones(size(Xs));
                    end    
            end
        end
        
        for j=1:length(Xs)
            if i == 1
                p2 = pdf('normal',Xs(j),a0,sqrt(sh));
            else
                Pa = pdf('normal',Xs(j),a1*Xs+a0,sqrt(sh));
                p2 = Pa* Pxz(i-1,:)';
            end
            Pxz(i,j)= p1(j) * p2;
        end
        Pxz(i,:) = Pxz(i,:)/sum(Pxz(i,:));
    end
    %% Smoothing section
%     Lxz=zeros(N,length(Xs));
%     for i=1:N-1
%         for j=1:length(Xs)
%            p2=0;
%            for k=1:length(Xs)
%                p2=p2+pdf('normal',Xs(j),a1*Xs(k)+a0,sqrt(sh))*Pxz(i,k);
%            end
%            Lxz(i+1,j)=p2;
%         end
%         Lxz(i+1,:)= Lxz(i+1,:)/sum(Lxz(i+1,:));
%     end
% 
%     Qxz= zeros(N,length(Xs));
%     Qxz(end,:)=Pxz(end,:);
%     for i=N-1:-1:1
%         for j=1:length(Xs)
%            p2=0;
%            for k=1:length(Xs)
%                p2=p2+pdf('normal',Xs(k),a1*Xs(j)+a0,sqrt(sh))*Qxz(i+1,k)/Lxz(i+1,k);
%            end
%            Qxz(i,j)=p2*Pxz(i,j);
%         end
%         Qxz(i,:) = Qxz(i,:)/sum(Qxz(i,:));
%     end
    %% Find Upper and Lower Bound for Both Filtering and Smoothing
%     for i=1:N
%         % mean
%         FMx(i)= Pxz(i,:)*Xs';
%         % Upper and Lower Bound for Filtering 
%         cs=cumsum(Pxz(i,:));
%         ind=find(cs>=0.025);
%         FLx(i)= Xs(ind(1)); 
%         ind=find(cs>=0.975);
%         FUx(i)= Xs(ind(1)); 
%         % mean
%         SMx(i)= Qxz(i,:)*Xs';
%         % Upper and Lower Bound for Filtering 
%         cs=cumsum(Qxz(i,:));
%         ind=find(cs>=0.025);
%         SLx(i)= Xs(ind(1)); 
%         ind=find(cs>=0.975);
%         SUx(i)= Xs(ind(1)); 
%     end
    
else % Gaussian approximate methods
    %% Run gaussian approximate for the Filtering part
    %----------------------------------------------%
    % XPrior
    XPre = zeros(N,1);
    % SPrior
    SPre = zeros(N,1);
    % FMxt
    FMx = zeros(N,1);
    % FSxt
    FSx = zeros(N,1);
    % Gamma type
    if OType==2
        for k=1:N
            % State update (As(1) is regressive term)
            if k == 1
                XPre(k) = a0;
                SPre(k) = sh;
            else
                XPre(k) = FMx(k-1)*a1 + a0;
                SPre(k) = a1*FSx(k-1)* a1 + sh;
            end
               
            if In(k)==1 % observed data
                % Observation: Gamma
                Mp  = exp( XPre(k)*b1 + b0);
                Yd  = Yn(k)- alpha;
                % Update
                 if ExistYb
                    Pk = exp(c2*XPre(k)+c1*exp(Yn(k))+c0)/(1+exp(c2*XPre(k)+c1*exp(Yn(k))+c0));
                    FSx(k)= 1/( (1/SPre(k)) + ((v*b1*b1*Yd)/Mp)+(c2*c2*(1-Pk)*Pk) );
                    FMx(k)= XPre(k) + FSx(k)* ( v * b1 *(Yd/Mp-1) + c2 * (Yb(k)-Pk));
                 else
                    FSx(k)= 1/((1/SPre(k)) + ((v*b1*b1*Yd)/Mp) );
                    FMx(k)= XPre(k) + FSx(k)* ( v * b1 *(Yd/Mp-1) );
                 end
            end
            if In(k)==0 && CType ==4
                % mean
                u = exp(b1*XPre(k)+b0);
                % probability of being above T
                pt    = gammainc(v*(T-alpha)/u,v,'upper');
                % pdf at T
                pf    = gamma_pdf(T-alpha,u,v);
                % de
                de    = pf - pt * v* b1*(((T-alpha)/u)-1);
                % update
                FSx(k)= 1/((1/SPre(k)) + (b1*(T-alpha) * pf * de / (pt*pt)) );
                FMx(k)= XPre(k) + FSx(k)* (b1 * (T-alpha) * pf/pt);
            end
%             if In(k)==0 && CType == 5    % Gamma
%                 min_x= XPre(k)- 100 * sqrt(SPre(k));
%                 max_x= XPre(k)+ 100 * sqrt(SPre(k));
%                 
%                 min_x= XPre(k)- 100 * sqrt(SPre(k));
%                 max_x= XPre(k)+ 100 * sqrt(SPre(k));
%                 %% 
%                 x  = min_x:0.01:max_x;
%                 u  = exp(b1.*x+b0);
%                 fx = gammainc(v*(T-alpha)./u,v,'upper').*exp(-(x-XPre(k)).*(x-XPre(k))/(2*SPre(k)));
%                 %c = integral(@(x)normal_c(x,b1,b0,so,XPre(k),SPre(k),T),min_x,max_x);
%                 %m = integral(@(x)normal_x(x,b1,b0,so,XPre(k),SPre(k),T),min_x,max_x);
%                 %s = integral(@(x)normal_xx(x,b1,b0,so,XPre(k),SPre(k),T),min_x,max_x);
%                 % update
%                 m  = (x*fx')/sum(fx);
%                 FMx(k)= m;
%                 FSx(k)= ((x-m).^2*fx')/sum(fx);
%                
%             end
        end
    end
    % Normal type
    if OType==1
        for k=1:N
            % State update (As(1) is regressive term)
            if k == 1
                XPre(k) = a0;
                SPre(k) = sh;
            else
                XPre(k) = FMx(k-1)*a1 + a0;
                SPre(k) = a1*FSx(k-1)* a1 + sh;
            end
                
            if In(k)==1 % observed data
                % Observation: Normal
                % Update
                 if ExistYb
                    Pk    = exp(c2*XPre(k)+c1*exp(Yn(k))+c0)/(1+exp(c2*XPre(k)+c1*exp(Yn(k))+c0));
                    FSx(k)= 1/((1/SPre(k)) + ((b1*b1)/so) +(c2*c2*(1-Pk)*Pk) );
                    FMx(k)= XPre(k) + FSx(k)* ( (b1*(Yn(k)-b1*(a1*XPre(k)+a0)-b0)/so) + c2*(Yb(k)-Pk));
                 else
                    FSx(k)= 1/((1/SPre(k)) + (b1*b1)/so );
                    FMx(k)= XPre(k) + FSx(k)* (b1*(Yn(k)-b1*(a1*XPre(k)+a0)-b0)/so);
                 end
            end
            if In(k)==0 && CType == 4
                % probability of being above T
                pt    = cdf('normal',T,b1*XPre(k)+b0,sqrt(so),'upper');
                % pdf at T
                pf    = pdf('normal',T,b1*XPre(k)+b0,sqrt(so));
                % de
                de    = pf - pt *(T-(b1*XPre(k)+b0))/so;
                FSx(k)= ((1/SPre(k)) + (b1*b1)*(pf/(pt*pt))*de )^-1;
                FMx(k)= XPre(k) + FSx(k)* (b1 * pf/pt);
            end
%             if In(k)==0 && CType == 5
%                 min_x= XPre(k)- 100 * sqrt(SPre(k));
%                 max_x= XPre(k)+ 100 * sqrt(SPre(k));
%                 %% 
%                 x  = min_x:0.01:max_x;
%                 fx = 0.5*(1-erf((T-b1*x-b0)/sqrt(2*so))).*exp(-(x-XPre(k)).*(x-XPre(k))/(2*SPre(k)));
%                 %c = integral(@(x)normal_c(x,b1,b0,so,XPre(k),SPre(k),T),min_x,max_x);
%                 %m = integral(@(x)normal_x(x,b1,b0,so,XPre(k),SPre(k),T),min_x,max_x);
%                 %s = integral(@(x)normal_xx(x,b1,b0,so,XPre(k),SPre(k),T),min_x,max_x);
%                 % update
%                 m  = (x*fx')/sum(fx);
%                 FMx(k)= m;
%                 FSx(k)= ((x-m).^2*fx')/sum(fx);
%             end
        end
    end
    
%     
%     %% Run the standard smoothing for the Smoothing part
%     %-------------------------------------%
%     % Kalman Smoothing
%     Ak        = zeros(N,1);
%     SMx      = zeros(N,1);
%     SMx(end) = FMx(end);
%     SSx      = zeros(N,1);
%     SSx(end) = FSx(end);
%     for k=N-1:-1:1
%         % Ak, equation (A.10)
%         Ak(k) = FSx(k) * a1 * ( SPre(k+1)^-1 );
%         % Smting function, equation (A.9)
%         SMx(k) = FMx(k) + Ak(k) * (SMx(k+1)- XPre(k+1));
%         % Variance update, equation (A.11)
%         SSx(k) = FSx(k) + Ak(k) * (SSx(k+1)- SPre(k+1)) * Ak(k);
%     end
%     
%     
%     Qxz = [];
%     Pxz = [];
%     FLx = FMx - 2 * sqrt(FSx);
%     FUx = FMx + 2 * sqrt(FSx);
%     SLx = SMx - 2 * sqrt(SSx);
%     SUx = SMx + 2 * sqrt(SSx);
%     

end


    %% nested function
    function f=gamma_pdf(y,u,v)
        f=((v*y/u)^v)*(1/y)*exp(-v*y/u)*(1/gamma(v));
    end
    %% Gamma function
    function f=gamma_c(x,b1,b0,v,alpha,xkk,sxx,T)
        u = exp(b1.*x+b0);
        f = gammainc(v*(T-alpha)./u,v,'upper').*exp(-(x-xkk).*(x-xkk)/(2*sxx))*(1/sqrt(2*pi*sxx));
    end
    function f=gamma_x(x,b1,b0,v,alpha,xkk,sxx,T)
        u = exp(b1.*x+b0);
        f = x.*gammainc(v*(T-alpha)./u,v,'upper').*exp(-(x-xkk).*(x-xkk)/(2*sxx))*(1/sqrt(2*pi*sxx));
    end
    function f=gamma_xx(x,b1,b0,v,alpha,xkk,sxx,T)
        u = exp(b1.*x+b0);
        f = x.*x.*gammainc(v*(T-alpha)./u,v,'upper').*exp(-(x-xkk).*(x-xkk)/(2*sxx))*(1/sqrt(2*pi*sxx));
    end
    function f=normal_c(x,b1,b0,so,xkk,sxx,T)
        f = 0.5*(1-erf((T-b1*x-b0)/sqrt(2*so))).*exp(-(x-xkk).*(x-xkk)/(2*sxx))*(1/sqrt(2*pi*sxx));
    end
    function f=normal_x(x,b1,b0,so,xkk,sxx,T)
        f = 0.5*x.*(1-erf((T-b1*x-b0)/sqrt(2*so))).*exp(-(x-xkk).*(x-xkk)/(2*sxx))*(1/sqrt(2*pi*sxx));
    end
    function f=normal_xx(x,b1,b0,so,xkk,sxx,T)
        f = 0.5*x.*x.*(1-erf((T-b1*x-b0)/sqrt(2*so))).*exp(-(x-xkk).*(x-xkk)/(2*sxx))*(1/sqrt(2*pi*sxx));
    end



    
end

