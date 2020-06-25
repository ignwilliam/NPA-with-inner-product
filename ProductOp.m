function [ Z ] = ProductOp( X,Y )
% takes product between operators X and Y, each of them are product of
% operators
% Z = X * Y

% define zero operator
zero.status = '0';
zero.as = '';
zero.ao = '';
zero.bs = '';
zero.bo = '';
zero.cs = '';
zero.co = '';

if strcmp(X.status,'0')||strcmp(Y.status,'0') % either X or Y is the zero operator
    Z = zero;
    return;
elseif strcmp(X.status,'I') % if X is the identity operator
    Z = Y;
    return;
elseif strcmp(Y.status,'I') % if Y is the identity operator
    Z = X;
    return;
else
    Z.status = '1';
    Z.as = strcat(X.as,Y.as);
    Z.ao = strcat(X.ao,Y.ao);
    Z.bs = strcat(X.bs,Y.bs);
    Z.bo = strcat(X.bo,Y.bo);
    Z.cs = strcat(X.cs,Y.cs);
    Z.co = strcat(X.co,Y.co);
    Z = CheckOrtho(Z); % check for orthogonality
    return;
end

