function v = WeightSort(k,r,g,y) %this function gives (1.4)
y=sort(y,'descend');
v1=y(1:k);
v2=y(k+1:length(y)).*(r/g);
v=[v1,v2];

end

