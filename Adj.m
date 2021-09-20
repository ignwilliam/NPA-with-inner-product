function [ Y ] = Adj( X )
% gives the adjoint of operators that are products of operators.
% input:
    % X: product of operators or an array of them
% output:
    % Y: adjoint of X

if numel(X) > 1
    % if X is an array of monomials
    dims = size(X);
    X = reshape(X,[1,numel(X)]);


    for i = 1:length(X)
        Y(i) = Adj(X(i));
    end

    Y = reshape(Y, dims);
else
    % the case if X is a single monomial
    if strcmp(X.status,'0') || strcmp(X.status,'I')
        % if X is a zero operator or identity operator
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

end

