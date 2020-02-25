function supp = Supp_calc(im,N)
    supp=zeros(N);
    for i=1:N
        for j=1:N
            if abs(im(i,j))>0
                supp(i,j)=1;
            end
        end
    end
end

