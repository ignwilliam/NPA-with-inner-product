function [ G, ref, I ] = GenerateMomentMatrix( S )
% generate moment matrix associated to a set of operators, S
% input:
    % S: a set of operators which are products of projectors
% output:
    % G: the moment matrix associated to set of operators, S

Id.status = 'I';
Id.as = []; Id.ao = [];
Id.bs = []; Id.bo = [];
Id.cs = []; Id.co = []; Id.cdagger = [];

ref = Id;


for i = 1:length(S)
    for j = 1:length(S)
        G(i,j) = ProductOp(Adj(S(i)),S(j));

        if strcmp(G(i,j).status, '0')
            I(i,j) = 0;
        else
            new_op = true;

            for k = 1:length(ref)
                if OpsCmp(G(i,j), ref(k))
                    I(i,j) = k;
                    new_op = false;
                    break;
                end
            end

            if new_op
                ref(end+1) = G(i,j);
                I(i,j) = length(ref);
            end
        end
    end
end


end

