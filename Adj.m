function [ Y ] = Adj( X )
% gives the adjoint of operators that are products of projectors.
% input:
    % X: product of projectors
% output:
    % Y: adjoint of X

if strcmp(X.status,'0') || strcmp(X.status,'I')
    Y = X;
    return;
else
    Y.status = '1';

    la = length(X.as);
    lb = length(X.bs);
    lc = length(X.cs);
    
    Y.as = flip(X.as);
    Y.ao = flip(X.ao);

    Y.bs = flip(X.bs);
    Y.bo = flip(X.bo);

    Y.cs = flip(X.cs);
    Y.co = flip(X.co);
    Y.cdagger = flip(X.cdagger);
    Y.cdagger = not(Y.cdagger);
end
end

