function [ G ] = GenerateMomentMatrix( S )
% generate moment matrix associated to a set of operators, S
% input:
    % S: a set of operators which are products of projectors
% output:
    % G: the moment matrix associated to set of operators, S

for i = 1:length(S)
    for j = 1:length(S)
        G(i,j) = ProductOp(Adj(S(i)),S(j));
    end
end


end

