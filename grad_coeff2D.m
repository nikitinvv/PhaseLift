function C = grad_coeff2D(X,N,M)
    Y=zeros(2*N-1);
    for p1=-N+1:N-1
        for p2=-N+1:N-1
            for q1=0:N-1
                for q2=0:N-1
                    
                    if 0 <= p1+q1 && p1+q1 <= N-1 && 0 <= p2+q2 && p2+q2 <= N-1
                        
                        Y(p1+N,p2+N)=Y(p1+N,p2+N)+X(q1+1+q2*N,q1+p1+1+(q2+p2)*N);
                        
                    end
                end
                
            end
        end
    end

%     Y2=zeros(2*N-1);
%     for p1=[-N+1:N-1]
%         for p2=[-N+1:N-1]
%             for q1=[0:N-1]
%                 for q2=[0:N-1]
%                     
%                     if 0 <= p1+q1 && p1+q1 <= N-1 && 0 <= p2+q2 && p2+q2 <= N-1
%                         
%                         Y2(p1+N,p2+N)=Y2(p1+N,p2+N)+X(MtV_ind([q1+1,q2+1],N),MtV_ind([q1+p1+1,q2+p2+1],N));
%                         
%                     end
%                 end
%                 
%             end
%         end
%     end
% 
% 
%     if(norm(Y-Y2)>0)
%         keyboard;
%     end    
        
    Y=padarray(Y,2*(M-N)+1,0,'post');
    Y=padarray(Y',2*(M-N)+1,0,'post')';
    
    C1=fft2(Y(1:M,1:M))+fft2(Y(M+1:2*M,1:M))+fft2(Y(1:M,M+1:2*M))+fft2(Y(M+1:2*M,M+1:2*M));
    
    
    C=zeros(N);
    for k1=[0:N-1]
        for k2=[0:N-1]
            C(k1+1,k2+1)=C1(k1+1,k2+1)*exp(-2*pi*i*(k2+k1)/M);
        end
    end
    C=conj(C);
    
    

end

