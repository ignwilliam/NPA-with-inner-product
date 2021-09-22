% test script for NPA with inner-product 
% security analysis for phase-encoded qubit BB84 via BFF bound

%% states and their Gram matrix
% phase-encoding BB84 states
% states = 1/sqrt(2) * [1, 1;
%                       1, -1;
%                       1, -1i;
%                       1, 1i];
states = [1, 0;
          0, 1;
          1/sqrt(2), 1/sqrt(2);
          1/sqrt(2), -1/sqrt(2)];
gram_input = GramMatrixQudit(states);

%% generate NPA relaxation
io_config = {[],[2,2],[2]};
level = 1;
extra_monomials = struct('status',{},'as',{},'ao',{},'bs',{},'bo',{},'cs',{},'co',{},'cdagger',{});
isComplex = false;
fprintf('Performing SDP relaxation...\n')
tic
npa = NPA(gram_input, io_config, level, extra_monomials, isComplex);
toc

%% BFF bound
num_nodes = 12;

%% compute phase-error rate
QBER_data = [];
R_data = [];
R_exact = [];

for err = 0:0.005:0.11
	tic
	H = BFF_bound(npa, err, num_nodes);
	toc
	R = H - h2(err);
	r = 1 - 2*h2(err);
	
	QBER_data = [QBER_data, err];	
	R_data = [R_data, R];
	R_exact = [R_exact, r];
end

plot(QBER_data,R_data, QBER_data, R_exact,'o')