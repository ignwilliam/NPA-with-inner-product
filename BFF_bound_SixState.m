function H = BFF_bound_SixState(npa, QBER, num_nodes)
	% compute lower bound on von Neumann entropy for phase-encoded qubit BB84
	% via BFF method and NPA with inner-product
	% inputs:
		% npa: the NPA relaxation object
		% QBER: the observed error rate (assuming depolarising noise)
		% num_nodes: number of nodes for the Gauss-Radau quadrature
	% output:
		% H: lower bound on the von Neumann entropy
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% formulate constraints
	% number of states
	num_state = length(npa.gram_input);

	% variables associated to Bob's moments
	B = GeneratePartyOps('B',3,1); % Bob's operators
	b = sdpvar(3,num_state,num_state);
	for i = 1:num_state
		for j = 1:num_state
			b(:,i,j) = npa.Operator2Variable(B,i,j);
		end
	end

	% %% QBER constraints:
	ez = QBER;
	ex = QBER;
	ey = QBER;

	% err_Z = 0.5 * ( 1 - b(1,1,1) + b(1,2,2) );
	% err_X = 0.5 * ( 1 - 2*real(b(2,1,2)) );
	% err_Y = 0.5 * ( 1 + 2*imag(b(3,1,2)) );
	% QBER_constraints = [err_Z <= ez, err_X <= ex, err_Y <= ey];


	%% full statistics constraints
	% for state |0>
	stat_constraints = [b(1,1,1) == 1-ez, b(2,1,1) == 1/2, b(3,1,1) == 1/2];

	% for state |1>
	stat_constraints = [stat_constraints, ...
						b(1,2,2) == ez, b(2,2,2) == 1/2, b(3,2,2) == 1/2];
	% for state |+>
	stat_constraints = [stat_constraints, ...
						0.5 * (1 + 2*real(b(1,2,1))) == 1/2, ...
						0.5 * (1 + 2*real(b(2,2,1))) == 1-ex, ...
						0.5 * (1 + 2*real(b(3,2,1))) == 1/2];
	% for state |->	
	stat_constraints = [stat_constraints, ...
						0.5 * (1 - 2*real(b(1,2,1))) == 1/2, ...
						0.5 * (1 - 2*real(b(2,2,1))) == ex, ...
						0.5 * (1 - 2*real(b(3,2,1))) == 1/2];
	% for state |+i>
	stat_constraints = [stat_constraints, ...
						0.5 * (1 - 2*imag(b(1,1,2))) == 1/2, ...
						0.5 * (1 - 2*imag(b(2,1,2))) == 1/2, ...
						0.5 * (1 - 2*imag(b(3,1,2))) == 1-ey];
	% for state |-i>	
	stat_constraints = [stat_constraints, ...
						0.5 * (1 + 2*imag(b(1,1,2))) == 1/2, ...
						0.5 * (1 + 2*imag(b(2,1,2))) == 1/2, ...
						0.5 * (1 + 2*imag(b(3,1,2))) == ey];
	
	% combine constraits
	constraints = [stat_constraints, npa.npa_constraints];
	
	%% BFF method
	% calculate the Gauss-Radau quadrature
	[t,w] = GenerateQuadrature(num_nodes);

	% calculate the lower bound on the von Neumann entropy
	H = 0; % total entropy
	for i = 1:num_nodes
		% calculate the quasi-relative entropy
		Hi = QuasiRelativeEntropy(npa,constraints,t(i),w(i));
		H = H + Hi;
	end

end