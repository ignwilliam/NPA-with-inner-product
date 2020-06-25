function [ Y ] = CheckOrtho( X )
% simplifies product between orthogonal projectors
% (e.g. A0 * A1 = 0, A0 * A0 = A0)
% input:
    % X: product of projectors
% output:
    % Y: simplified operator
    
zero.status = '0'; %define zero
zero.as = '';
zero.ao = '';
zero.bs = '';
zero.bo = '';
zero.cs = '';
zero.co = '';

dummy = X; %as X changes during the process, dummy is fixed reference of X

if ~isempty(dummy.as)
    for k = length(dummy.as)-1:-1:1
        if strcmp(X.as(k),X.as(k+1)) %check for same inputs of Alice
            if ~strcmp(X.ao(k),X.ao(k+1)) %check for different outputs
                Y = zero;
                return;
            else
               X.as(k+1) = '';
               X.ao(k+1) = '';
            end
        end
    end
else
end

if ~isempty(dummy.bs)
    for k = length(dummy.bs)-1:-1:1 %Same for Bob.
        if strcmp(X.bs(k),X.bs(k+1))
           if ~strcmp(X.bo(k),X.bo(k+1))
              Y = zero;
              return;
           else
               X.bs(k+1) = '';
               X.bo(k+1) = '';
           end
        end
    end
end

if ~isempty(dummy.cs)
    for k = length(dummy.cs)-1:-1:1 %Same for Charlie
        if strcmp(X.cs(k),X.cs(k+1))
           if ~strcmp(X.co(k),X.co(k+1))
              Y = zero;
              return;
           else
               X.cs(k+1) = '';
               X.co(k+1) = '';
           end
        end
    end
end

Y = X;

end

