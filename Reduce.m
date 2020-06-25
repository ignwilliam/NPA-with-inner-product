function X = Reduce(Y) 
% takes a set of operators and remove the repeated elements

% removing zeros
for i = length(Y):-1:1
    if strcmp(Y(i).status, '0')
        Y(i) = [];
    end
end

%removing repeated operators
for i = 1:length(Y)
    for j = length(Y):-1:i+1
        if isequal(Y(i), Y(j))
            Y(j) = [];
        end
    end
end
X = Y;