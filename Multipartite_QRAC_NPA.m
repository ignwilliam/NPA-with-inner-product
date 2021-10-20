% QRAC with two receivers (based on arXiv: 1803.04796)

% we use the NPA method with inner product constraints to bound the winning
% probability of Alice when Bob's winning probability is fixed

% the referee uses the BB84 states:
    % |z = 00> = |0>
    % |z = 01> = |->
    % |z = 10> = |+>
    % |z = 11> = |1>

% the winning probabilities are defined as P(a = z_x) and P(b = z_y)
% where x and y are inputs of Alice and Bob while
% a and b are their outputs

% IMPORTANT:    
% here, we use CVX for disciplined convex programming with MOSEK
% solver refer to http://cvxr.com/cvx/doc/mosek.html for an instruction on 
% how to use MOSEK with CVX

%% Bob's winning probability
N = 20;
t = linspace(0,pi/2,N);
R = (cos(pi/8)^2 - 0.5);
pb = 0.5 + R * sin(t);
% pb = linspace(0.5,1/2 * (1 + 1/sqrt(2)),N); % Bob's winning probability
pa = zeros(1,N);

%% define the P&M network scenario
nX = 2; % number of inputs for Alice
nA = 2; % number of outputs for Alice
nY = 2; % number of inputs for Bob
nB = 2; % number of outputs for Bob
nZ = 0; % number of inputs for Charlie
nC = 0; % number of outputs for Charlie

%% the level of the NPA hierarchy
Q = 1;

%% generate moment matrix
S = GenerateOps(nX,nA,nY,nB,nZ,nC,Q); % the operators for level-Q NPA hierarchy
ref = GenerateOps(nX,nA,nY,nB,nZ,nC,2*Q); % reference table
G = GenerateMomentMatrix(S); % NPA moment matrix
I = GenerateIndexMatrix(G,ref); % give indices to the elements of G

L = length(S);

%% relevant operators
% identity operators
Id.status = 'I';
Id.as = '';  Id.ao = '';
Id.bs = '';  Id.bo = '';
Id.cs = '';  Id.co = '';

% Alice's projectors
A00.status = '1';
A00.as = '0'; A00.ao = '0';
A00.bs = '';  A00.bo = '';
A00.cs = '';  A00.co = '';

A10.status = '1';
A10.as = '1'; A10.ao = '0';
A10.bs = '';  A10.bo = '';
A10.cs = '';  A10.co = '';

% Bob's projectors
B00.status = '1';
B00.as = '';  B00.ao = '';
B00.bs = '0'; B00.bo = '0';
B00.cs = '';  B00.co = '';

B10.status = '1';
B10.as = '';  B10.ao = '';
B10.bs = '1'; B10.bo = '0';
B10.cs = '';  B10.co = '';

% AB operators
A00B00.status = '1';
A00B00.as = '0'; A00B00.ao = '0';
A00B00.bs = '0'; A00B00.bo = '0';
A00B00.cs = '';  A00B00.co = '';

A10B00.status = '1';
A10B00.as = '1'; A10B00.ao = '0';
A10B00.bs = '0'; A10B00.bo = '0';
A10B00.cs = '';  A10B00.co = '';

A00B10.status = '1';
A00B10.as = '0'; A00B10.ao = '0';
A00B10.bs = '1'; A00B10.bo = '0';
A00B10.cs = '';  A00B10.co = '';

A10B10.status = '1';
A10B10.as = '1'; A10B10.ao = '0';
A10B10.bs = '1'; A10B10.bo = '0';
A10B10.cs = '';  A10B10.co = '';

%% Gram matrix of the input states
lambda = [         1,  1/sqrt(2), 1/sqrt(2),          0;
           1/sqrt(2),          1,         0, -1/sqrt(2);
           1/sqrt(2),          0,         1,  1/sqrt(2);
                   0, -1/sqrt(2), 1/sqrt(2),          1];


%% solve SDP
for n = 1:N
cvx_begin sdp quiet
    cvx_solver mosek
    variable v(length(ref),4,4);
    expression pA; % Alice's winning probability
    expression pB; % Bob's winning probability
    expression g(4,4,L,L); % moment matrices   

    % projectors
    expression id;
    expression a00;
    expression a10;
    expression b00;
    expression b10;
    expression a00b00;
    expression a10b00;
    expression a00b10;
    expression a10b10;

    for i = 1:4
        for j = 1:4
            % assign variable to projectors
            id(i,j)  = v(IndexOps(Id,ref),i,j);
            a00(i,j) = v(IndexOps(A00,ref),i,j);
            a10(i,j) = v(IndexOps(A10,ref),i,j);
            b00(i,j) = v(IndexOps(B00,ref),i,j);
            b10(i,j) = v(IndexOps(B10,ref),i,j);

            a00b00(i,j) = v(IndexOps(A00B00,ref),i,j);
            a10b00(i,j) = v(IndexOps(A10B00,ref),i,j);
            a00b10(i,j) = v(IndexOps(A00B10,ref),i,j);
            a10b10(i,j) = v(IndexOps(A10B10,ref),i,j);
        end
    end

    % assign variable to moment matrix
    for i = 1:4
        for j = 1:4
            for k = 1:L
                for l = 1:L
                    if I(k,l) == 0
                        g(k,l,i,j) = 0;
                    else
                        g(k,l,i,j) = v(I(k,l),i,j);
                    end
                end
            end
        end
    end
    
    % the big Gram matrix
    gamma = [ g(:,:,1,1), g(:,:,1,2), g(:,:,1,3), g(:,:,1,4);
              g(:,:,2,1), g(:,:,2,2), g(:,:,2,3), g(:,:,2,4);
              g(:,:,3,1), g(:,:,3,2), g(:,:,3,3), g(:,:,3,4);
              g(:,:,4,1), g(:,:,4,2), g(:,:,4,3), g(:,:,4,4)];

    % define winning probabilities
    pA = 1/8 * ( 4 + a00(1,1) + a10(1,1) + a00(2,2) - a10(2,2) - a00(3,3) + a10(3,3) - a00(4,4) - a10(4,4) );
    pB = 1/8 * ( 4 + b00(1,1) + b10(1,1) + b00(2,2) - b10(2,2) - b00(3,3) + b10(3,3) - b00(4,4) - b10(4,4) );

    maximize pA;
    subject to
        % overlap constraint
        for i = 1:4
            for j = 1:4
                id(i,j) == lambda(i,j);
            end
        end

        % PSD constraint
        gamma >= 0;

        % Bob's winning probability constraint
        pB == pb(n);
cvx_end

pa(n) = pA;
end

x = linspace(0.5,1/2 * (1 + 1/sqrt(2)),1000);
y = 0.5 + 0.5 * sqrt( 0.5 - (2 *x - 1).^2);

plot(x,y,pb,pa,'o');