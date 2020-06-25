function [ z ] = OpsCmp( X , Y )
% compares 2 operators X and Y and give 1 if they are the same and 0
% if they are not

% based on our code, two operators are equivalent iff all the fields of the
% structures that represent them are the same

%% compare each field of the structure
z1 = strcmp(X.status,Y.status);
z2 = strcmp(X.as,Y.as);
z3 = strcmp(X.ao,Y.ao);
z4 = strcmp(X.bs,Y.bs);
z5 = strcmp(X.bo,Y.bo);
z6 = strcmp(X.cs,Y.cs);
z7 = strcmp(X.co,Y.co);

% if all fields are the same, then the sum will be 7
if z1 + z2 + z3 + z4 + z5 + z6 + z7 == 7
    z = 1;
    return;
else
    z = 0;
    return;
end



end

