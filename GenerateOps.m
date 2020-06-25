function ops = GenerateOps(nX,nA,nY,nB,nZ,nC,Q)

% generate set of operators for Alice, Bob, Charlie to build the moment
% matrix for a given level of the hierarchy Q

% inputs:
    % nX: number of inputs for Alice
    % nA: number of outputs for Alice
    % nY: number of inputs for Bob
    % nB: number of outputs for Bob
    % nZ: number of inputs for Charlie
    % nC: number of outputs for Charlie
    % Q: level of the hierarchy
% output:
    % ops: set of operators associated to the moment matrix

% NOTES:
% the building blocks are projective measurements
% we use structures to represent the operators, the data type for each
% field is string that represents the inputs and outputs of each party (from left to right)
    % status: indicates whether the operator is zero operator, identity operator or product of projectors
    % as: the settings of Alice's operators
    % ao: the outputs of Alice's operators
    % bs: the settings of Bob's operators
    % bo: the outputs of Bob's operators   
    % cs: the settings of Charlie's operators
    % co: the outputs of Charlie's operators

% example 1:
    % A00 is Alice's projector corresponding to input 0 and ouput 0
    % then we have:
    % A00.status = '1'
    % A00.as = '0'
    % A00.ao = '0'
    % A00.bs = ''
    % A00.bo = ''
    % A00.cs = ''
    % A00.co = ''
    
% example 2:
    % A00B10 is an operator which corresponds to product of Alice's and
    % Bob's projectors when Alice's input is 0 and her ouput is 0,
    % while Bob's input is 1 and his output is 0,
    % then we have:
    % A00B10.status = '1'
    % A00B10.as = '0'
    % A00B10.ao = '0'
    % A00B10.bs = '1'
    % A00B10.bo = '0'
    % A00B10.cs = ''
    % A00B10.co = ''

% example 3:
    % A00A10 is an operator which corresponds to product of two of Alice's
    % projectors: A00 * A10 (the first projector corresponds to when her
    % input is 0 and output is 0, while the second projector is when her
    % input is 1 and output is 1)
    % then we have:
    % A00A10.status = '1'
    % A00A10.as = '01'
    % A00A10.ao = '00'
    % A00A10.bs = ''
    % A00A10.bo = ''
    % A00A10.cs = ''
    % A00A10.co = ''

S = []; % initialise the set

%% the identity operator
S(1).status = 'I'; % status: '0' is zero, '1' is non-empty, 'I' is identity
S(1).as = ''; % Alice's setting
S(1).ao = ''; % Alice's output
S(1).bs = ''; % Bob's setting
S(1).bo = ''; % Bob's output
S(1).cs = ''; % Charlie's setting
S(1).co = ''; % Charlie's output

P = []; % Q1 operators
i = 1;

%% Alice's operators
for x = 1:nX
    for a = 1:nA-1
        P(i).status = '1'; % status: '0' is zero, '1' is non-empty, 'I' is identity
        P(i).as = int2str(x-1); % Alice's setting
        P(i).ao = int2str(a-1); % Alice's output
        P(i).bs = ''; % Bob's setting
        P(i).bo = ''; % Bob's output
        P(i).cs = ''; % Charlie's setting
        P(i).co = ''; % Charlie's output
        i = i+1;
    end
end

%% Bob's operators
for y = 1:nY
    for b = 1:nB-1
        P(i).status = '1'; % status: '0' is zero, '1' is non-empty, 'I' is identity
        P(i).as = ''; % Alice's setting
        P(i).ao = ''; % Alice's output
        P(i).bs = int2str(y-1); % Bob's setting
        P(i).bo = int2str(b-1); % Bob's output
        P(i).cs = ''; % Charlie's setting
        P(i).co = ''; % Charlie's output
        i = i+1;
    end
end

%% Charlie's operators
for z = 1:nZ
    for c = 1:nC-1
        P(i).status = '1'; % status: '0' is zero, '1' is non-empty, 'I' is identity
        P(i).as = ''; % Alice's setting
        P(i).ao = ''; % Alice's output
        P(i).bs = ''; % Bob's setting
        P(i).bo = ''; % Bob's output
        P(i).cs = int2str(z-1); % Charlie's setting
        P(i).co = int2str(c-1); % Charlie's output
        i = i+1;
    end
end

%% operators for higher moments
S = [S,P];

for l = 2:Q
    L = length(S);
    for j = 1:length(P)
        for k = 1:L
            S(length(S)+1) = ProductOp(S(k),P(j)); 
        end
    end
    S = Reduce(S); % remove identical operators from the set S
end

ops = S;

end

