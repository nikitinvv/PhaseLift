function X=Astarnomasks2d(c,n)%Computes A*. Takes a matrix of size mxm where m geq n
m=length(c(:,1));
chat=fft2(c);
%chat(2:end,2:end)=chat(2:end,2:end);



%this loop creates the toeplitz matrices that will go in the larger block
%toeplitz matrix X
for j=1:m
rows=chat(1:n,j);
%we want the last element on the first upper diagonal after main diagonal.
%This is the last number in chat, which explains the fliplr. chat(1) is
%added since otherwise toeplitz complains of diagonal conflict.
columns=chat(m-n+2:m,j);columns=[columns ; chat(1,j)];columns=flipud(columns);
Y(j)={toeplitz(rows,columns)};
end
%we now do an index matrix so that the Ys end up in the correct position,
%the idea behind the indexing is as in line 12
ind=toeplitz([1,m:-1:m-n+2],[1:n])';
X=cell2mat(Y(ind));