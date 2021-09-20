function ops = GenerateOps(io_config, level)

% generate set of operators for Alice, Bob, Eve to build the moment
% matrix for a given level of the hierarchy
% we assume that Alice's and Bob's operators are projective measurements
% and Eve's operators are arbitrary operators


% inputs:
    % io_config: input-output configuration
    %            (a cell array that specifies the number of outputs for each input for each party)
    % level: relaxation level of the hierarchy
% output:
    % ops: set of operators associated to the moment matrix

% NOTES:
% we use structures to represent the operators, the data type for each
% field is an array that represents the inputs and outputs of each party (from left to right)
    % status: indicates whether the operator is zero operator ('0'), identity operator ('I') or product of operators ('1')
    % as: the settings of Alice's operators
    % ao: the outputs of Alice's operators
    % bs: the settings of Bob's operators
    % bo: the outputs of Bob's operators   
    % cs: the settings of Eve's operators
    % co: the outputs of Eve's operators
    % cdagger: boolean variable that indicates whether the operator is hermitian conjugated

% example 1:
    % A00 is Alice's projector corresponding to input 0 and ouput 0
    % then we have:
    % A00.status = '1'
    % A00.as = 0
    % A00.ao = 0
    % A00.bs = []
    % A00.bo = []
    % A00.cs = []
    % A00.co = []
    % A00.cdagger = []
    
% example 2:
    % A00B10 is an operator which corresponds to product of Alice's and
    % Bob's projectors when Alice's input is 0 and her ouput is 0,
    % while Bob's input is 1 and his output is 0,
    % then we have:
    % A00B10.status = '1'
    % A00B10.as = 0
    % A00B10.ao = 0
    % A00B10.bs = 1
    % A00B10.bo = 0
    % A00B10.cs = []
    % A00B10.co = []

% example 3:
    % A00A10 is an operator which corresponds to product of two of Alice's
    % projectors: A00 * A10 (the first projector corresponds to when her
    % input is 0 and output is 0, while the second projector is when her
    % input is 1 and output is 1)
    % then we have:
    % A00A10.status = '1'
    % A00A10.as = [0,1]
    % A00A10.ao = [0,0]
    % A00A10.bs = []
    % A00A10.bo = []
    % A00A10.cs = []
    % A00A10.co = []

% number of settings
nX = length(io_config{1}); nY = length(io_config{2}); nZ = length(io_config{3});


S = []; % initialise the set

%% the identity operator
S(1).status = 'I'; % status: '0' is zero, '1' is non-empty, 'I' is identity
S(1).as = []; % Alice's setting
S(1).ao = []; % Alice's output
S(1).bs = []; % Bob's setting
S(1).bo = []; % Bob's output
S(1).cs = []; % Eve's setting
S(1).co = []; % Eve's output
S(1).cdagger = []; % Eve's dagger

P = []; % Q1 operators
i = 1;

%% Alice's operators
for x = 1:nX
    nA = io_config{1}(x);
    for a = 1:nA-1
        P(i).status = '1'; % status: '0' is zero, '1' is non-empty, 'I' is identity
        P(i).as = x-1; % Alice's setting
        P(i).ao = a-1; % Alice's output
        P(i).bs = []; % Bob's setting
        P(i).bo = []; % Bob's output
        P(i).cs = []; % Eve's setting
        P(i).co = []; % Eve's output
        P(i).cdagger = [];
        i = i+1;
    end
end

%% Bob's operators
for y = 1:nY
    nB = io_config{2}(y);
    for b = 1:nB-1
        P(i).status = '1'; % status: '0' is zero, '1' is non-empty, 'I' is identity
        P(i).as = []; % Alice's setting
        P(i).ao = []; % Alice's output
        P(i).bs = y-1; % Bob's setting
        P(i).bo = b-1; % Bob's output
        P(i).cs = []; % Eve's setting
        P(i).co = []; % Eve's output
        P(i).cdagger = [];
        i = i+1;
    end
end

%% Eve's operators
for z = 1:nZ
    nC = io_config{3}(z);
    for c = 1:nC
        for dagger = 0:1
            P(i).status = '1'; % status: '0' is zero, '1' is non-empty, 'I' is identity
            P(i).as = []; % Alice's setting
            P(i).ao = []; % Alice's output
            P(i).bs = []; % Bob's setting
            P(i).bo = []; % Bob's output
            P(i).cs = z-1; % Eve's setting
            P(i).co = c-1; % Eve's output
            P(i).cdagger = dagger;
            i = i+1;
        end
    end
end

%% operators for higher moments
S = [S,P];

for l = 2:level
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

