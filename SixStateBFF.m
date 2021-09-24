% test script for NPA with inner-product 
% security analysis for phase-encoded qubit six-state protocol via BFF bound

%% states and their Gram matrix
states = [1, 0;
          0, 1];
gram_input = GramMatrixQudit(states);

%% generate NPA relaxation
io_config = {[],[2,2,2],[2]};
level = 1;
extra_monomials = struct('status',{},'as',{},'ao',{},'bs',{},'bo',{},'cs',{},'co',{},'cdagger',{});
isComplex = true;
fprintf('Performing SDP relaxation...\n')
tic
npa = NPA(gram_input, io_config, level, extra_monomials, isComplex);
toc

%% BFF bound
num_nodes = 8;

% save data
QBER_data = [];
R_data = [];
R_exact = [];

% collecting data
tic
for err = 0:0.006:0.126
	H = BFF_bound_SixState(npa, err, num_nodes); % calculate von Neumann entropy
	R = H - h2(err);
	
	QBER_data = [QBER_data, err];	
	R_data = [R_data, R];
end
toc

QBER = linspace(0,0.126,1E4);
R_BB84 = 1 - 2.*h2(QBER);
R_SS = 1 - h2(QBER) - QBER - (1-QBER).*h2( (1-3*QBER/2) ./ (1-QBER));

plot(QBER_data,R_data, 'o', QBER, R_BB84, QBER, R_SS)