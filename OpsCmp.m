function [ z ] = OpsCmp( X , Y )
% compares operators X and Y and give 1 if they are the same and 0
% if they are not

% based on our code, two operators are equivalent iff all the fields of the
% structures that represent them are the same

%% compare each field of the structure
eq_array = [];

eq_array(end+1) = strcmp(X.status, Y.status);
eq_array(end+1) = isequal(X.as, Y.as);
eq_array(end+1) = isequal(X.ao, Y.ao);
eq_array(end+1) = isequal(X.bs, Y.bs);
eq_array(end+1) = isequal(X.bo, Y.bo);
eq_array(end+1) = isequal(X.cs, Y.cs);
eq_array(end+1) = isequal(X.co, Y.co);
eq_array(end+1) = isequal(X.cdagger, Y.cdagger);

z = all(eq_array);

end

