function iotaPhaseLift2D(M,N,N_masks,n_trials,noise_min,noise_step,noise_max,Nest,x_0M)

%HERE M IS THE NUMBER OF ROWS, CHANGE OF MEANING WRT 1D CASE
opts = optimset('Diagnostics','off', 'Display','off');

noise=[];
noise1=[];
rank_stat1=[];

for w=[noise_min:noise_step:noise_max]
    trials=[];
    trials1=[];
    for tr=1:n_trials
        %x_0M=randn(M,N)+i*randn(M,N);
        
        %x_0M=x_0M'; %<-
        x_0=x_0M(:);
        x_I=x_0;
        x_0=x_0/norm(x_0);
        X_0=x_0*x_0';
        
        
        
        dist_from_ground=[];
        eigen=[];
        subsequent=[];
        
        lent=size(X_0);
        len=lent(1);
        meas_data=[];
        
        storeMasks=cell(1,N_masks+1);
        storediagMasks=cell(1,N_masks+1);
        for v=2:N_masks+1
            storeMasks{v}=randi([0 1],M,N);
            storediagMasks{v}=diag(storeMasks{v}(:));
        end
        storeMasks{1}=ones(N);
        storediagMasks{1}=eye(N^2);
        
        meas_data=[];
        for l=1:N_masks+1
            meas_data=[meas_data;abs(grad_coeff2D(storediagMasks{l}*X_0*storediagMasks{l},N,N))];
        end
        norm(meas_data)
        
        temp_noise=randn(N*(N_masks+1),N);
        temp_noise=w*norm(meas_data,'fro')*temp_noise/norm(temp_noise,'fro');
        meas_data=real(meas_data)+temp_noise; %kill the imaginary part that MATLAB might keep track of; meas_data must be a real vector
        
        
        
        
        
        
        t=1/(N_masks*M^3.5+1); %stepsize definition
        gamma=N^4;
        
        X_k1temp=zeros(N^2);
        X0=eye(N^2);
        X_k=X_k1temp;
        X_k1=X_k1temp;
        X_k0=X0;
        r=0;
        
        while (r<20000 && norm(X_k0-X_k1temp)>1e-10)
            
            grad1=0;
            r=r+1;
            
            %Nesterov acceleration
            if Nest
                
                X_k1=X_k1temp + (r-1)/(r+2)*(X_k1temp - X_k0);
                X_k0=X_k1temp;
            else
                X_k1=X_k1temp;
            end
            
            c=cell(1,N_masks+1);
            d=cell(1,N_masks+1);
            
            for p=1:N_masks+1
                c{p}=zeros(N);
                c{p}=(grad_coeff2D(storediagMasks{p}*X_k1*storediagMasks{p},N,N) - meas_data(N*(p-1)+1:N*p,:));
            end
            
            
            
            for b=1:N_masks+1
                grad1=grad1+GRAD2(c{b},N,storeMasks{b});
            end
            
            
            
            Xhat1 = X_k1 - t*grad1;
            
            %             if r>1 && norm(prox1(2:end))<1e-15
            %                 [P1,D1]=eigs(Xhat1,1,'LM');
            %                 X_k1temp=P1*D1*P1';
            %
            %             else
            
            [P1,D1]=eig(Xhat1);
            prox1=(real(diag(D1))/t-gamma*Prox1(1,gamma,t,real(diag(D1))')')/(1/t - gamma);
            X_k1temp=P1*diag(prox1)*P1';
            %             end
            
            
            if ismember(r,[200:200:20000])
                [x_k1,ek1] = eigs(X_k1temp,1);
                fun1=@(c)minERROR(c,x_0,x_k1);
                phase1=fsolve(fun1,[0;0;1],opts);
                B=norm(x_I)*reshape(x_k1/(phase1(1)+1i*phase1(2)),[N,N]);
                figure
                imshow(B,'InitialMagnification',1000)
                
                
                norm(B-x_0M)
                
                dist_from_ground=[dist_from_ground,norm(X_k1temp-X_0)]
                t_eig=eig(X_k1temp);
                eigen=[eigen,t_eig(1:10)]
                subsequent=[subsequent,norm(X_k0-X_k1temp)]
                
                
            end
            norm(X_0-X_k1temp)
        end
        
        
        
        % fun=@(c)minERROR(c,x_0,x_k);
        % phase=fsolve(fun,[0;0;1],opts);
        [x_k1,ek1] = eigs(X_k1temp,1);
        fun1=@(c)minERROR(c,x_0,x_k1);
        phase1=fsolve(fun1,[0;0;1],opts);
        norm((phase1(1)+1i*phase1(2))*x_0-x_k1)^2/norm(x_0)^2;
        
        % norm_t=[norm((phase(1)+1i*phase(2))*x_0-x_k)^2,norm((phase1(1)+1i*phase1(2))*x_0-x_k1)^2,norm(x_k*x_k'-x_0*x_0','fro'),norm(x_k1*x_k1'-x_0*x_0','fro')]
        % error_track=[error_track,norm_t];
        
        
        
        A=reshape(x_k1/(phase1(1)+1i*phase1(2)),[N,N]);
        figure(1)
        imshow(norm(x_I)*A)
        figure(2)
        imshow(x_0M)
        norm(norm(x_I)*A-x_0M)
        
        
        %vrank_stat1=[rank_stat1,rank(X_k1,1e-3)]
    end
    
    
    
    noise=[noise,sum(trials)/n_trials]
    noise1=[noise1,sum(trials1)/n_trials]
end

