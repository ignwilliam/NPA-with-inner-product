% Compute BB84 with BFF method

% we use the NPA method with inner product constraints to bound the 
% von Neumann entropy for BB84 protocol

% IMPORTANT:    
% here, we use CVX for disciplined convex programming with MOSEK
% solver refer to http://cvxr.com/cvx/doc/mosek.html for an instruction on 
% how to use MOSEK with CVX


%% define the P&M network scenario
nX = 0; % number of inputs for Alice
nA = 0; % number of outputs for Alice
nY = 2; % number of inputs for Bob
nB = 2; % number of outputs for Bob
nZ = 1; % number of inputs for Charlie
nC = 2; % number of outputs for Charlie

%% the level of the NPA hierarchy
Q = 2;

%% generate moment matrix
S = GenerateOps(nX,nA,nY,nB,nZ,nC,Q); % the operators for level-Q NPA hierarchy
ref = GenerateOps(nX,nA,nY,nB,nZ,nC,2*Q); % reference table
G = GenerateMomentMatrix(S); % NPA moment matrix
I = GenerateIndexMatrix(G,ref); % give indices to the elements of G

L = length(S);

% %% relevant operators
% identity operator
id.status = 'I';
id.as = []; id.ao = [];
id.bs = []; id.bo = [];
id.cs = []; id.co = []; id.cdagger = [];

% Bob's operators
B00.status = '1';
B00.as = []; B00.ao = [];
B00.bs = 0;  B00.bo = 0;
B00.cs = []; B00.co = []; B00.cdagger = [];


B10.status = '1';
B10.as = []; B10.ao = [];
B10.bs = 1;  B10.bo = 0;
B10.cs = []; B10.co = []; B10.cdagger = [];

% Z operators
Z0.status = '1';
Z0.as = []; Z0.ao = [];
Z0.bs = []; Z0.bo = [];
Z0.cs = 0;  Z0.co = 0; Z0.cdagger = 0;

Z0d.status = '1';
Z0d.as = []; Z0d.ao = [];
Z0d.bs = []; Z0d.bo = [];
Z0d.cs = 0;  Z0d.co = 0; Z0d.cdagger = 1;

Z1.status = '1';
Z1.as = []; Z1.ao = [];
Z1.bs = []; Z1.bo = [];
Z1.cs = 0;  Z1.co = 1; Z1.cdagger = 0;

Z1d.status = '1';
Z1d.as = []; Z1d.ao = [];
Z1d.bs = []; Z1d.bo = [];
Z1d.cs = 0;  Z1d.co = 1; Z1d.cdagger = 1;

% ZZ term
Z0Z0d = ProductOp(Z0,Z0d);
Z0dZ0 = ProductOp(Z0d,Z0);
Z1Z1d = ProductOp(Z1,Z1d);
Z1dZ1 = ProductOp(Z1d,Z1);

%% Gram matrix of the input states
states = [1, 0;
          0, 1;
          1/sqrt(2), 1/sqrt(2);
          1/sqrt(2), -1/sqrt(2)];

lambda = GramMatrixQudit(states);


% %% solve SDP
% cvx_begin sdp
%     cvx_solver mosek
%     variable v(length(ref),4,4);
%     expression pA; % Alice's winning probability
%     expression pB; % Bob's winning probability
%     expression g(4,4,L,L); % moment matrices   

%     % projectors
%     expression id;
%     expression a00;
%     expression a10;
%     expression b00;
%     expression b10;
%     expression a00b00;
%     expression a10b00;
%     expression a00b10;
%     expression a10b10;

%     for i = 1:4
%         for j = 1:4
%             % assign variable to projectors
%             id(i,j)  = v(IndexOps(Id,ref),i,j);
%             a00(i,j) = v(IndexOps(A00,ref),i,j);
%             a10(i,j) = v(IndexOps(A10,ref),i,j);
%             b00(i,j) = v(IndexOps(B00,ref),i,j);
%             b10(i,j) = v(IndexOps(B10,ref),i,j);

%             a00b00(i,j) = v(IndexOps(A00B00,ref),i,j);
%             a10b00(i,j) = v(IndexOps(A10B00,ref),i,j);
%             a00b10(i,j) = v(IndexOps(A00B10,ref),i,j);
%             a10b10(i,j) = v(IndexOps(A10B10,ref),i,j);
%         end
%     end

%     % assign variable to moment matrix
%     for i = 1:4
%         for j = 1:4
%             for k = 1:L
%                 for l = 1:L
%                     if I(k,l) == 0
%                         g(k,l,i,j) = 0;
%                     else
%                         g(k,l,i,j) = v(I(k,l),i,j);
%                     end
%                 end
%             end
%         end
%     end
    
%     % the big Gram matrix
%     gamma = [ g(:,:,1,1), g(:,:,1,2), g(:,:,1,3), g(:,:,1,4);
%               g(:,:,2,1), g(:,:,2,2), g(:,:,2,3), g(:,:,2,4);
%               g(:,:,3,1), g(:,:,3,2), g(:,:,3,3), g(:,:,3,4);
%               g(:,:,4,1), g(:,:,4,2), g(:,:,4,3), g(:,:,4,4)];

%     % define winning probabilities
%     pA = 1/8 * ( 4 + a00(1,1) + a10(1,1) + a00(2,2) - a10(2,2) - a00(3,3) + a10(3,3) - a00(4,4) - a10(4,4) );
%     pB = 1/8 * ( 4 + b00(1,1) + b10(1,1) + b00(2,2) - b10(2,2) - b00(3,3) + b10(3,3) - b00(4,4) - b10(4,4) );

%     maximize pA;
%     subject to
%         % overlap constraint
%         for i = 1:4
%             for j = 1:4
%                 id(i,j) == lambda(i,j);
%             end
%         end

%         % PSD constraint
%         gamma >= 0;

%         % Bob's winning probability constraint
%         pB == pb;
% cvx_end