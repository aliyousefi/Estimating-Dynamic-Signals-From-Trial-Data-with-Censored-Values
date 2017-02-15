function [a,b] = ay_censored_data()
% censored data filtering
close all
% Numebr of samples
N  = 100;

% Hidden state model (we assume x0 is 0)
a   = 0.5;
sx  = 0.2025;
Xn  = zeros(N,1);
Xn(1) = sqrt(sx)*randn;
for i=2:N
    Xn(i)=Xn(i-1)*a + sqrt(sx)*randn;
end
% Observation 
Zn  = zeros(N,1);
b   = 0.4;
sy  = 0.001;
Zn  = b * Xn + sqrt(sy) * ones(size(Xn));

% Observation 

%% The Main Filtering Algorithm
xs  = -1.0:0.05:1.0;
T   = 0.8;
for t = 1:4
    T   = T / 2;
    % Find missing points
    In  = ones(N,1);
    ind = find(Zn>T);
    In(ind)=0;
    % Calculate the Exact posterior
    Pxz= zeros(N,length(xs));
    for i=1:N
        for j=1:length(xs)
            if In(i)==1 % normal
                p1 = pdf('normal',Zn(i),b*xs(j),sqrt(sy));
            else
                p1 = 1-cdf('normal',T,b*xs(j),sqrt(sy));
            end
            p2 = 0;
            if i == 1
                    p2 = pdf('normal',xs(j),0,sqrt(sx));
            else
                for k=1:length(xs)
                    p2 = p2 + Pxz(i-1,k)*pdf('normal',xs(j),a*xs(k),sqrt(sx));
                end
            end

            Pxz(i,j)= p1 * p2;
        end
        Pxz(i,:) = Pxz(i,:)/sum(Pxz(i,:));
        
    end
    
    % Calculate approximate posterior
    Mx = zeros(N,1);
    Sx = zeros(N,1);
    for k=1:N
        if k == 1
            XPre = 0;
            SPre = sx;
        else
            XPre = a * Mx(k-1);
            SPre = a*Sx(k-1)*a+ sx;
        end
        a1 = pdf('normal',T,b*XPre,sqrt(sy));
        a2 = cdf('normal',T,b*XPre,sqrt(sy),'upper');
        a3 = (T-b*XPre)/sy;
        Sx(k) = (1/SPre + In(k)*(b*b/sy) + (1-In(k))* ((b*b*a1)/(a2*a2))*(a1-a2*a3))^-1;
        Mx(k) = XPre + Sx(k) * ( (In(k)*b * (Zn(k)-b*XPre)/sy) + ((1-In(k))*b*a1/a2) );
    end
    
    figure(t);
    subplot(2,1,1);
    %imagesc(log(Pxz'))
    imagesc(Pxz')
    text=['T=' num2str(T)];
    title(text);
    xlabel('K');
    set(gca,'YTick',linspace(1,length(xs),7));
    %set(gca,'YTickLabel',[-1.5 -1.0 -0.5 0 0.5 1.0 1.5])
    set(gca,'YTickLabel',[-1.0 -0.66 -0.33 0 0.33 0.66 1])
    hold on;
    ind=find(In==0);
    plot(ind,In(ind),'y*');
    hold off;
    
    subplot(2,1,2);
    Mn = zeros(N,1);
    Ss = zeros(N,1);
    for i=1:N
        Mn(i) = sum(Pxz(i,:).*xs);
        Ss(i) = sqrt(sum((xs-Mn(i)).*(xs-Mn(i)).*Pxz(i,:)));
    end
    ay_plot_bound(1, 1:N ,Mn, (Mn-2*Ss)', (Mn+2*Ss)');
    hold on
    %ay_plot_bound(1, 1:N ,Mx, (Mx-2*sqrt(Sx))', (Mx-2*sqrt(Sx))');
    plot(1:N,Mx,'g','LineWidth',2);
    hold on;
    plot(Xn,'r');
    ind=find(In==0);
    plot(ind,In(ind),'r*');
    hold off;
    text=['T=' num2str(T)];
    title(text);
    xlabel('K');
 
    figure(4+t)
    surf(Pxz');
    text=['T=' num2str(T)];
    title(text);
    set(gca,'YTick',linspace(1,length(xs),7));
    set(gca,'YTickLabel',[-1.0 -0.66 -0.33 0 0.33 0.66 1])
    xlabel('n')
    ylabel('X')
    axis([1 N 1  length(xs) 0 Inf])
end

end

