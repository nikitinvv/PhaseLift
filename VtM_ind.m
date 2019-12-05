function p_m = VtM_ind(p,N)
    p1=mod(p,N)+(mod(p,N)==0);
    p2=floor(p/N)+(~(mod(p,N)==0));
    p_m=[p1,p2];
end

