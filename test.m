% test script for NPA with inner-product 
% security analysis for phase-encoded qubit BB84 via phase-error rate bound

%% set observed QBER
QBER = 0.01;

%% states and their Gram matrix
% phase-encoding BB84 states
states = 1/sqrt(2) * [1, 1;
                      1, -1;
                      1, -1i;
                      1, 1i];

gram_input = GramMatrixQudit(states);

%% generate NPA relaxation
io_config = {[],[2,2],[]};
level = 1;
extra_monomials = struct('status',{},'as',{},'ao',{},'bs',{},'bo',{},'cs',{},'co',{},'cdagger',{});
isComplex = true;
npa = NPA(gram_input, io_config, level, extra_monomials, isComplex);

%% compute phase-error rate
QBER_data = [];
eph_data = [];
R_data = [];
R_analytical = [];

for err = 0:0.005:0.11
	ep = PhaseError(npa,err,gram_input);
	R = 1 - h2(ep) - h2(err);
	K = 1 - 2*h2(err);
	QBER_data = [QBER_data, err];
	eph_data = [eph_data, ep];
	R_data = [R_data, R];
	R_analytical = [R_analytical, K];
end

plot(QBER_data,R_data, QBER_data, R_analytical,'o')