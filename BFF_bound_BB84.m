function H = BFF_bound_BB84(npa, QBER, num_nodes)
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
	B = GeneratePartyOps('B',2,1); % Bob's operators
	b = sdpvar(2,num_state,num_state);
	for i = 1:num_state
		for j = 1:num_state
			b(:,i,j) = npa.Operator2Variable(B,i,j);
		end
	end

	%% QBER constraints:
	% % for (X,Y) basis protocol
	% ex = QBER;
	% ey = QBER;
	% err_X = 0.5 * ( 1-b(1,1,1) + b(1,2,2) );
	% err_Y = 0.5 * ( 1-b(2,3,3) + b(2,4,4) );
	% QBER_constraints = [err_X <= ex, err_Y <= ey];

	% for (Z,X)-basis protocol
	ez = QBER;
	ex = QBER;
	err_Z = 0.5 * ( 1 - b(1,1,1) + b(1,2,2) );
	err_X = 0.5 * ( 1 - 2*real(b(2,1,2)) );
	QBER_constraints = [err_Z <= ez, err_X <= ex];

	% combine constraits
	constraints = [QBER_constraints, npa.npa_constraints];
	
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