%The function inserts columns and rows of zeros, in the positions given by
%zero_pos
function [padded_X] = pseudo_padding(X,zero_pos)
    padded_X=X;
        for i=zero_pos
            padded_X(i,:)=0;
            padded_X(:,i)=0;
        end
        


end

