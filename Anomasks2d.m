function c=Anomasks2d(X,n,m)%this function computes A(X), or, which is the same, the coefficients c in the paper
%X should be an n^2 matrix, m as in paper
Y=zeros(2*m,2*m);%Y has to be bigger than 2n for the modulo operations in the end
X=reshape(X,[n,n,n,n]);%turns the matrix into a 4 tensor as in the paper, we can now work with the multiindices as in the paper, 
%except for shifts of +1 and +m+1 in different places...
for p1=-n+1:n-1
    for p2=-n+1:n-1
        P1=p1+m+1;P2=p2+m+1;%for indexing and zero padding purposes the vector Y is indexed form P=1:2m, instead of p=-n+1:n-1 as in the paper.
        for q1=max(-p1,0):min(n-1,n-1-p1)
            for q2=max(-p2,0):min(n-1,n-1-p2)
                Y(P1,P2)=Y(P1,P2)+X(q1+1,q2+1,p1+q1+1,p2+q2+1);%this is the q summation in formula (17)
            end
        end
    end
end
Ymod=Y(1:m,1:m)+Y(1:m,m+1:2*m)+Y(m+1:2*m,1:m)+Y(m+1:2*m,m+1:2*m);%p index 0 corresponds to P=m+1. 
%This is the modulo operation performed before (18)
c=conj(fft2(Ymod));
