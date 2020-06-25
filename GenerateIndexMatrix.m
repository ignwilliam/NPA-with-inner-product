function I = GenerateIndexMatrix(G,ref)
% generate a matrix of indices of the elements of the moment matrix which
% inputs:
    % G: moment matrix
    % ref: a reference table to identify each element of the moment matrix
% output:
    % I: matrix of indices

d = size(G,1);
I = zeros(d);

for i = 1:d
    for j = 1:d
        I(i,j) = IndexOps(G(i,j),ref); % assign index to each element of G
    end
end