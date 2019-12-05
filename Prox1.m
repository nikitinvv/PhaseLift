function m = Prox1(k,g,t,y) %k is the cardinality, r=rho and g=gamma, y the point in which we compute prox.
[y1,I]=sort(y,'descend');

Z=zeros(length(y));
for h=1:length(y)
    Z(h,I(h))=1;
end

if  y1(k)<0    %if the vector is already "suitable" I return 1.4; note that in 1.4 k still depends on y!
    d=find(y1<0,1);
    v1=y1(1:d);
    v2=y1(d+1:length(y1))*(1/(g*t));
    m=(Z'*[v1 v2]')';
    
    
elseif (1/(g*t))*y1(k+1)<y1(k)
    m=(Z'*sort(WeightSort(k,1/t,g,y1),'descend')')';
    
else %We compute the indices called l* and j* in the paper; we call them just l and j.
    y=WeightSort(k,1/t,g,y1);
    j=1;
    while y(j)>y(k+1)
        j=j+1;
    end
    
    l=k;
    while  y(l+1)>=y(k)
        l=l+1;
        if l==length(y)
            break
        end
    end
    
    z=sort(y(j:l),'descend'); %We compute the vector z according to point 3.
    for i=1:length(z)-1
        s=(z(i)+z(i+1))/2;
        x=Candidate(k,s,y);
        j1=1; %we compute the indices called l and j in the paper. We call them l1 and j1.
        while x(j1)~=s
            j1=j1+1;
        end
        l1=length(y);
        while x(l1)~=s
            l1=l1-1;
        end
        
        sI=((1/t)*sum(y1(j1:l1)))/((k+1-j1)*(1/t) + (l1-k)*g);
        if (sI>=z(i+1) && z(i)>=sI)
            m=(Z'*sort(Candidate(k,sI,y),'descend')')';
            
            break
        end
    end
    
    
end



