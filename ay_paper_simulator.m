function ay_paper_simulator()

T= linspace(0.1,1.5,2);
Xs = -4:0.2:4;

fid = fopen('normal.txt','wt');
fclose(fid);

fid = fopen('normal_count.txt','wt');
fclose(fid);


fid = fopen('normal_binary.txt','wt');
fclose(fid);

fid = fopen('normal_binary_count.txt','wt');
fclose(fid);
    
for rpt = 1:10
    % create data
    [Xk,Yn,Yb,In,Param] =ay_data_generator(100,1,1,0.95,0.08,T);
    
    for t=1:length(T)
        text=['Iteration = ' num2str(rpt) ', Threshold= '  num2str(T(t))]; 
        disp(text)
        
        pnt = length(In(:,t))-sum(In(:,t));
        % now run different methods(exact, just normal)
        Pxz = ay_smoothing(1,1,Yn,[],In(:,t),Xs,Param,T(t));
        mse1 = ay_mse(Xk,Xs,Pxz);
        cnt1 = ay_count(Xk,Xs,Pxz);
        


        % now run different methods(impute, just normal)
        Pxz =ay_smoothing(1,2,Yn,[],In(:,t),Xs,Param,T(t));
        mse2 = ay_mse(Xk,Xs,Pxz);
        cnt2 = ay_count(Xk,Xs,Pxz);


        % now run different methods(ignor, just normal)
        Pxz =ay_smoothing(1,3,Yn,[],In(:,t),Xs,Param,T(t));
        mse3 = ay_mse(Xk,Xs,Pxz);
        cnt3 = ay_count(Xk,Xs,Pxz);

        % now run different methods(Gaussian approximate, just normal)
        [temp,Mx,Sx] =ay_smoothing(1,4,Yn,[],In(:,t),Xs,Param,T(t));
        mse4 = mean((Xk-Mx).^2);
        cnt4 = ay_count(Xk,Mx,Sx,1);
        
        fid = fopen('normal.txt','at');
        fprintf(fid,'%f  %f  %f  %f  ',mse1,mse2,mse3,mse4);
        fclose(fid);
    
        fid = fopen('normal_count.txt','at');
        fprintf(fid,'%f  %f  %f  %f  %f  ',cnt1,cnt2,cnt3,cnt4,pnt);
        fclose(fid);
    

        % now run different methods(exact, just normal)
        Pxz = ay_smoothing(1,1,Yn,[],In(:,t),Xs,Param,T(t));
        mse1 = ay_mse(Xk,Xs,Pxz);
        cnt1 = ay_count(Xk,Xs,Pxz);

        % now run different methods(impute, just normal)
        Pxz =ay_smoothing(1,2,Yn,[],In(:,t),Xs,Param,T(t));
        mse2 = ay_mse(Xk,Xs,Pxz);
        cnt2 = ay_count(Xk,Xs,Pxz);

        % now run different methods(ignor, just normal)
        Pxz =ay_smoothing(1,3,Yn,[],In(:,t),Xs,Param,T(t));
        mse3 = ay_mse(Xk,Xs,Pxz);
        cnt3 = ay_count(Xk,Xs,Pxz);

        % now run different methods(Gaussian approximate, just normal)
        [temp,Mx,Sx] =ay_smoothing(1,4,Yn,Yb,In(:,t),Xs,Param,T(t));
        mse4 = mean((Xk-Mx).^2);
        cnt4 = ay_count(Xk,Mx,Sx,1);
        
        
        fid = fopen('normal_binary.txt','at');
        fprintf(fid,'%f  %f  %f  %f  ',mse1,mse2,mse3,mse4);
        fclose(fid);
    
        fid = fopen('normal_count_binary.txt','at');
        fprintf(fid,'%f  %f  %f  %f  %f  ',cnt1,cnt2,cnt3,cnt4,pnt);
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
