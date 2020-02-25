N=64;
nu=1;

opts = optimset('Diagnostics','off', 'Display','off');

%here we generate "sparse" ground truths. This will encode the
%knowledge of the support, i.e. we know where the measured object
%is not

%random ground truth and then zero-padded and lifted
M=imread('Mango64.jpg');
M=im2double(M);

x_0=M(:,:,2);
%x_0=padarray(x_0,[N/2 N/2],0,'both');
pdd_im=x_0;



%support definition
Z=Supp_calc(x_0,N);
true_supp=Z;
Z=Z(:);
g_const=Z*Z';

g_truth=x_0;
n_g_truth=norm(g_truth,'fro');

x_0=x_0(:);
norm_x_0=norm(x_0);
x_0=x_0/norm(x_0);
X_0=x_0*x_0';



meas_data=real(Anomasks2d(X_0,N,2*N));
b=norm(meas_data,'fro');

gamma=N;
t=1/N^5;

%starting points
X0=zeros(N^2);
X_k=X0;
Z_k=X0;
Y_k=X0;
for e=1:5
    tic
    [P1,D1]=eig(Z_k-Y_k*t);
    prox1=(real(diag(D1))/t-gamma*Prox1(1,gamma,t,real(diag(D1))')')/(1/t - gamma);
    
    X_k=P1*diag(prox1)*P1';
    
    %FBS routine to compute Z_k
    Z_k1=zeros(N^2);
    Z_k0=Z_k1;
    r=1;
    tau=1/((N^2)^2+1/t);
    while norm(Z_k1-Z_k0,'fro')>1e-2 || r<2
        
        Z_k1=Z_k1 + (r-1)/(r+2)*(Z_k1 - Z_k0);
        Z_k0=Z_k1;
        r=r+1;
        Z_hat=Z_k1-tau*((Z_k1-X_k-Y_k*t)/t + Astarnomasks2d(Anomasks2d(real(Z_k1),N,2*N)-meas_data,N));
        Z_k1=Z_hat.*g_const;
        
        
        %                 Z_k0=Z_k1;
        %                 Z_k1=Z_hat;
        
    end
    Z_k=Z_k1;
    %norm(Anomasks2d(Z_k1,N,2*N)-meas_data,'fro')/b
    
    Y_k=Y_k+(X_k-Z_k)/t;
    [x_k1,ek1] = eigs(X_k,1);
    
    norm(Z_k1-X_0,'fro');
    rsh=padarray(reshape(x_k1,N,N),[N/2 N/2],0,'both');
    norm(meas_data-abs(fft2(rsh)).^2,'fro');
    toc
    
    %                         fun=@(c)minERROR(c,x_0,x_k1);
    %                         phaseAft1=fsolve(fun,[0;0;0],opts)
    %norm(x_0-abs(x_k1))
    
    
    %imshow(norm_x_0*rec,'InitialMagnification',1000)
    
    
    
    
end

[x_k1,ek1] = eigs(Z_k,1);











exit()
