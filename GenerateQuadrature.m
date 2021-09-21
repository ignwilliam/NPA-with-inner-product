
function [t,w] = GenerateQuadrature(n)
% Function to compute nodes & weights of Radau quadrature
% with uniform weight in the interval [0,1] (the fixed node is 1).

% Written by Wen Yu Kon
% inputs:
    % n: number of nodes for Gauss-Radau quadrature
% outputs:
    % t: the nodes of the Radau quadrature
    % w: the corresponding weight of each node

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % nitialise nodes & weights (For [-1,1] case).
    t1 = zeros(1,n-1);
    w1 = zeros(1,n-1);

    % Initialise to find some points of the function.
    points = 100*n;
    val = zeros(1,points);

    % Find the value of these points.
    x_val = linspace(-1+1/points,1-1/points,points);
    for i = 1:points
        val(i) = fun(n,x_val(i));
    end

    % Find out where the points change in sign.
    val = sign(val);
    val_change = diff(val);

    % Indices where values change (i/i+1).
    ind = find(val_change ~= 0);

    % Boundaries to seek.
    x_lb = x_val(ind);
    x_ub = x_val(ind+1);

    % Flag if there is errors.
    if length(x_lb) ~= n-1
        t = [];
        w = [];
        fprintf('Error in Radau points!\n')
    end

    % Function.
    fun2 = @(x) fun(n,x);

    % Compute nodes and weights.
    for i = 1:n-1
        t1(i) = fzero(fun2,[x_lb(i),x_ub(i)]);
        w1(i) = 1/n^2 * (1-t1(i)) / (leg_poly(n-1,t1(i)))^2;
    end

    % Adjustment to [0,1] case.
    t = 0.5 * (1-t1);
    w = 0.5 * w1;

    % Include the fixed node
    t = [1, t];
    w = [1/n^2, w];

    % Sort nodes in ascending order
    t = flip(t);
    w = flip(w);
end

%Function used to find the nodes.
function out = fun(n,x)
val = legendre(n,x);
val2 = legendre(n-1,x);
out = (val(1) + val2(1))/(1+x);
end

%Legendre polynomial.
function out = leg_poly(n,x)
val = legendre(n,x);
out = val(1);
end