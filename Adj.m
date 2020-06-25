function [ Y ] = Adj( X )
% gives the adjoint of operators that are products of projectors.
% input:
    % X: product of projectors
% output:
    % Y: adjoint of X

if strcmp(X.status,'0')||strcmp(X.status,'I')
    Y = X;
    return;
else
    Y.status = '1';

    la = length(X.as);
    lb = length(X.bs);
    lc = length(X.cs);
    
    % reverse the order for Alice
    if la>1
        for k=1:la
            Y.as(k) = X.as(la+1-k);
            Y.ao(k) = X.ao(la+1-k);
        end
    else
        Y.as = X.as;
        Y.ao = X.ao;
    end
    % reverse the order for Bob
    if lb>1
        for k=1:lb
            Y.bs(k) = X.bs(lb+1-k);
            Y.bo(k) = X.bo(lb+1-k);
        end
    else
        Y.bs = X.bs;
        Y.bo = X.bo;
    end
    
    % reverse the order for Charlie
    if lc>1
        for k=1:lc
            Y.cs(k) = X.cs(lc+1-k);
            Y.co(k) = X.co(lc+1-k);
        end
    else
        Y.cs = X.cs;
        Y.co = X.co;
    end

end
end

