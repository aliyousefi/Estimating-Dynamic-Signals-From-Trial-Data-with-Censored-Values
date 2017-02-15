function ay_paper_simulator_parallel()

T= linspace(log(0.4),log(2.4),18);
Xs = -2.5:0.02:3.5;

fid = fopen('normal.txt','wt');
fclose(fid);

fid = fopen('normal_count.txt','wt');
fclose(fid);


fid = fopen('normal_binary.txt','wt');
fclose(fid);

fid = fopen('normal_binary_count.txt','wt');
fclose(fid);
    
for rpt = 1:500
    % create data
    [Xk,Yn,Yb,In,Param] =ay_data_generator(100,1,0.25,0.95,0.02,T);
    
    parfor t=1:length(T)
        text=['Iteration = ' num2str(rpt) ', Threshold= '  num2str(T(t))]; 
        disp(text)
        
        pnt(t) = sum(In(:,t));
        % now run different methods(exact, just normal)
        
        Pxz = ay_smoothing(1,1,Yn,[],In(:,t),Xs,Param,T(t));
        mse1(t) = ay_mse(Xk,Xs,Pxz);
        cnt1(t) = ay_count(Xk,Xs,Pxz);
        

        % now run different methods(impute, just normal)
        %disp('IMPUTE ')
        Pxz =ay_smoothing(1,2,Yn,[],In(:,t),Xs,Param,T(t));
        mse2(t) = ay_mse(Xk,Xs,Pxz);
        cnt2(t) = ay_count(Xk,Xs,Pxz);


        % now run different methods(ignore, just normal)
        Pxz =ay_smoothing(1,3,Yn,[],In(:,t),Xs,Param,T(t));
        mse3(t) = ay_mse(Xk,Xs,Pxz);
        cnt3(t) = ay_count(Xk,Xs,Pxz);

        % now run different methods(Gaussian approximate, just normal)
        [temp,Mx,Sx] =ay_smoothing(1,4,Yn,[],In(:,t),Xs,Param,T(t));
        for s=1:100
            Pxz(s,:)=exp(-(Xs-Mx(s)).^2/(2*Sx(s)));
            Pxz(s,:)=Pxz(s,:)/sum(Pxz(s,:));
        end
        mse4(t) = ay_mse(Xk,Xs,Pxz);
        cnt4(t) = ay_count(Xk,Xs,Pxz);
        
        % now run different methods(impute, just normal)
        %disp('MULTIPLE IMPUTE ')
        Pxz =ay_smoothing(1,5,Yn,[],In(:,t),Xs,Param,T(t));
        mse5(t) = ay_mse(Xk,Xs,Pxz);
        cnt5(t) = ay_count(Xk,Xs,Pxz);

        
     
        % now run different methods(exact, just normal)
        Pxz = ay_smoothing(1,1,Yn,Yb,In(:,t),Xs,Param,T(t));
        mse6(t) = ay_mse(Xk,Xs,Pxz);
        cnt6(t) = ay_count(Xk,Xs,Pxz);

        % now run different methods(impute, just normal)
        %disp('IMPUTE  MIX')
        Pxz =ay_smoothing(1,2,Yn,Yb,In(:,t),Xs,Param,T(t));
        mse7(t) = ay_mse(Xk,Xs,Pxz);
        cnt7(t) = ay_count(Xk,Xs,Pxz);

        % now run different methods(ignor, just normal)
        Pxz =ay_smoothing(1,3,Yn,Yb,In(:,t),Xs,Param,T(t));
        mse8(t) = ay_mse(Xk,Xs,Pxz);
        cnt8(t) = ay_count(Xk,Xs,Pxz);

        % now run different methods(Gaussian approximate, just normal)
        [temp,Mx,Sx] =ay_smoothing(1,4,Yn,Yb,In(:,t),Xs,Param,T(t));
        for s=1:100
            Pxz(s,:)=exp(-(Xs-Mx(s)).^2/(2*Sx(s)));
            Pxz(s,:)=Pxz(s,:)/sum(Pxz(s,:));
        end
        mse9(t) = ay_mse(Xk,Xs,Pxz);
        cnt9(t) = ay_count(Xk,Xs,Pxz);
        
        % now run different methods(impute, just normal)
        %disp('MULTIPLE IMPUTE  MIX')
        Pxz =ay_smoothing(1,5,Yn,Yb,In(:,t),Xs,Param,T(t));
        mse10(t) = ay_mse(Xk,Xs,Pxz);
        cnt10(t) = ay_count(Xk,Xs,Pxz);

        
    end
    
    for t=1:length(T)
        fid = fopen('normal.txt','at');
        fprintf(fid,'%f  %f  %f  %f  %f ',mse1(t),mse2(t),mse3(t),mse4(t),mse5(t));
        fclose(fid);
    
        fid = fopen('normal_count.txt','at');
        fprintf(fid,'%f  %f  %f  %f  %f  %f ',cnt1(t),cnt2(t),cnt3(t),cnt4(t),cnt5(t),pnt(t));
        fclose(fid);
    
        fid = fopen('normal_binary.txt','at');
        fprintf(fid,'%f  %f  %f  %f  %f ',mse6(t),mse7(t),mse8(t),mse9(t),mse10(t));
        fclose(fid);

        fid = fopen('normal_binary_count.txt','at');
        fprintf(fid,'%f  %f  %f  %f  %f  %f ',cnt6(t),cnt7(t),cnt8(t),cnt9(t),cnt10(t),pnt(t));
        fclose(fid);
    end
    
    fid = fopen('normal.txt','at');
    fprintf(fid,'\r\n');
    fclose(fid);
    
    fid = fopen('normal_count.txt','at');
    fprintf(fid,'\r\n');
    fclose(fid);
    
        
    fid = fopen('normal_binary.txt','at');
    fprintf(fid,'\r\n');
    fclose(fid);
    
    fid = fopen('normal_binary_count.txt','at');
    fprintf(fid,'\r\n');
    fclose(fid);
    
end    
