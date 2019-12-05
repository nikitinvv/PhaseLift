function w = Candidate(k,s,y) %this function gives (1.6)
w=zeros(1,length(y));
for i=1:k
    w(i)=max(y(i),s);
end
for i=k+1:length(y)
    w(i)=min(y(i),s);
end
end

