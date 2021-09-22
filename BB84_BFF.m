% test script for NPA with inner-product 
% security analysis for phase-encoded qubit BB84 via BFF bound

%% states and their Gram matrix
% phase-encoding BB84 states
% states = 1/sqrt(2) * [1, 1;
%                       1, -1;
%                       1, -1i;
%                       1, 1i];
states = [1, 0;
          0, 1];
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
num_nodes = 8;

% save data
QBER_data = [];
R_data = [];
R_exact = [];

% collecting data
tic
for err = 0:0.005:0.11
	H = BFF_bound_BB84(npa, err, num_nodes); % calculate von Neumann entropy
	R = H - h2(err);
	
	QBER_data = [QBER_data, err];	
	R_data = [R_data, R];
end
toc

QBER = linspace(0,0.11,1E4);
R_exact = 1 - 2.*h2(QBER);
plot(QBER_data,R_data, 'o', QBER, R_exact)