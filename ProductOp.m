function [ Z ] = ProductOp( X,Y )
% takes product between operators X and Y, each of them are product of
% operators
% Z = X * Y

% define zero operator
zero.status = '0';
zero.as = [];
zero.ao = [];
zero.bs = [];
zero.bo = [];
zero.cs = [];
zero.co = [];
zero.cdagger = [];

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
    Z.as = [X.as,Y.as];
    Z.ao = [X.ao,Y.ao];
    Z.bs = [X.bs,Y.bs];
    Z.bo = [X.bo,Y.bo];
    Z.cs = [X.cs,Y.cs];
    Z.co = [X.co,Y.co];
    Z.cdagger = [X.cdagger,Y.cdagger];
    Z = CheckOrtho(Z); % check for orthogonality
    return;
end

