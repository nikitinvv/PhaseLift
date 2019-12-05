function F = minERROR(c,x0,xk)
    len=size(xk);
    xkD=zeros(len(1),1);
    x0D=zeros(len(1),1);
    
        for k=1:len(1)
            xkD(k)=(c(1)*real(x0(k))-c(2)*imag(x0(k))-real(xk(k)))*real(x0(k)) + (c(1)*imag(x0(k))+c(2)*real(x0(k))-imag(xk(k)))*imag(x0(k));
            x0D(k)=(c(1)*real(x0(k))-c(2)*imag(x0(k))-real(xk(k)))*(-imag(x0(k))) + (c(1)*imag(x0(k))+c(2)*real(x0(k))-imag(xk(k)))*real(x0(k));
        end
        
        F(1)=sum(xkD)+c(1)*c(3);
        F(2)=sum(x0D)+c(2)*c(3);
        F(3)=c(1)^2+c(2)^2-1;
end

