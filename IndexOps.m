function [ Y ] = IndexOps( Op, Ref )
% identifies an operator by searching through a reference table for the
% same operator

% inputs:
    % Op: the operator to be identified
    % Ref: the reference table of all operators in the moment matrix
% output:
    % Y: an index number that identifies the operator

if strcmp(Op.status,'0')
    Y = 0;
else
    for i = 1:length(Ref)
        if OpsCmp(Op,Ref(i))
            Y = i;
            return
        end
    end
end

end
