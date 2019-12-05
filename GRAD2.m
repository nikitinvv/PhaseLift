function final_grad2 = GRAD2(C,N,storeM)
    T=fft2(C);
    Ma=find(storeM);
    final_grad2=zeros(N^2); %I write N here, but this has to be changed if we want oversampling!!!
    for im =1:length(Ma)
        [m1,m2]=ind2sub([N,N],Ma(im));
        m1 = m1-1;
        m2 = m2-1;        
        for in =1:length(Ma)
            [n1,n2]=ind2sub([N,N],Ma(in));
            n1 = n1-1;
            n2 = n2-1;
            if m1-n1 >= 0 && m2-n2 >= 0
                final_grad2(Ma(im),Ma(in))=T(m1-n1+1,m2-n2+1);
            elseif m1-n1 < 0 && m2-n2 < 0
                final_grad2(Ma(im),Ma(in))=T(N+m1-n1+1,N+m2-n2+1);
            elseif m1-n1 < 0 && m2-n2 >= 0
                final_grad2(Ma(im),Ma(in))=T(N+m1-n1+1,m2-n2+1);
            elseif m1-n1 >= 0 && m2-n2 < 0
                final_grad2(Ma(im),Ma(in))=T(m1-n1+1,N+m2-n2+1);
            end           
        end
    end
    
%     final_grad=zeros(N^2); %I write N here, but this has to be changed if we want oversampling!!!
%         for m1=[0:N-1] 
%             for m2=[0:N-1]
%                 for n1=[0:N-1]
%                     for n2=[0:N-1]
%                          if ismember(MtV_ind([m1+1,m2+1],N),Ma) && ismember(MtV_ind([n1+1,n2+1],N),Ma)
%                             
%                             if m1-n1 >= 0 && m2-n2 >= 0
%                                 final_grad(MtV_ind([m1+1,m2+1],N),MtV_ind([n1+1,n2+1],N))=T(m1-n1+1,m2-n2+1);
%                             elseif m1-n1 < 0 && m2-n2 < 0
%                                 final_grad(MtV_ind([m1+1,m2+1],N),MtV_ind([n1+1,n2+1],N))=T(N+m1-n1+1,N+m2-n2+1);
%                             elseif m1-n1 < 0 && m2-n2 >= 0
%                                 final_grad(MtV_ind([m1+1,m2+1],N),MtV_ind([n1+1,n2+1],N))=T(N+m1-n1+1,m2-n2+1);
%                             elseif m1-n1 >= 0 && m2-n2 < 0
%                                 final_grad(MtV_ind([m1+1,m2+1],N),MtV_ind([n1+1,n2+1],N))=T(m1-n1+1,N+m2-n2+1);
%                             end
%                         end
%                      end
%                 end
%             end
%         end
% 
% % check         
%     if(norm(final_grad-final_grad2)>0)
%       keyboard;    
%     end
    
end

