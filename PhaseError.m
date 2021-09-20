function [ep] = PhaseError(npa, QBER, gram_input)
	% calculate the phase-error rate for 
	% qubit phase-encoding BB84 protocol via SDP relaxation

	% inputs:
		% npa: the NPA relaxation with inner-product
		% QBER: the error rate (we assume depolarising channel)
		% Gram input: the Gram matrix of the input states
	% output:
		% ep: phase-error rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% number of states
	num_state = length(gram_input);

	% variables associated to Bob's moments
	B = GeneratePartyOps('B',2,1); % Bob's operators
	b = sdpvar(2,num_state,num_state);
	for i = 1:num_state
		for j = 1:num_state
			b(:,i,j) = npa.Operator2Variable(B,i,j);
		end
	end

	% objective function:
	eph = 0.5 - imag(b(2,1,2)); % phase error rate

	% QBER constraints:
	ex = QBER;
	ey = QBER;
	err_X = 0.5 * ( 1-b(1,1,1) + b(1,2,2) );
	err_Y = 0.5 * ( 1-b(2,3,3) + b(2,4,4) );
	QBER_constraints = [err_X <= ex, err_Y <= ey];

	% combine constraits
	constraints = [QBER_constraints, eph <= 0.5, npa.npa_constraints];

	num_constraints = length(constraints)

	options = sdpsettings('solver', 'mosek', 'verbose', 0);
	sol = optimize(constraints, -eph, options);
	ep = value(eph);
end