function Hi = QuasiRelativeEntropy(npa,constraints,ti,wi)
	% compute the quasi-relative entropies using BFF bound
	% inputs:
		% npa: the NPA relaxation object
		% constraints: the constraints for the optimisation
		% ti: the i-th node of the Gauss-Radau quadrature
		% wi: the weight of the i-th node
	% output:
		% Hi: the quasi-relative entropy

	%% find variables that form the objective function
	% number of states
	num_state = length(npa.gram_input);

	% relevant operators in computing the quasi-relative entropy
	B = GeneratePartyOps('B',2,1); % Bob's operators
	Z = GeneratePartyOps('Z',1,2)'; % Eve's operators
	Zd = Adj(Z); % the adjoint

	% product of Eve's operators
	for i = 1:2
		ZZd(i) = ProductOp(Z(i), Zd(i));
		ZdZ(i) = ProductOp(Zd(i), Z(i));
	end

	% the associated variables
	z = sdpvar(2,num_state,num_state);
	zd = sdpvar(2,num_state,num_state);
	zzd = sdpvar(2,num_state,num_state);
	zdz = sdpvar(2,num_state,num_state);

	for i = 1:num_state
		for j = 1:num_state
			z(:,i,j) = npa.Operator2Variable(Z,i,j);
			zd(:,i,j) = npa.Operator2Variable(Zd,i,j);
			zzd(:,i,j) = npa.Operator2Variable(ZZd,i,j);
			zdz(:,i,j) = npa.Operator2Variable(ZdZ,i,j);
		end
	end

	%% formulate objective function and solve SDP
	obj = 0.5 * real(z(1,1,1) + zd(1,1,1) + (1-ti)*zdz(1,1,1) + ...
			     	 z(2,2,2) + zd(2,2,2) + (1-ti)*zdz(2,2,2) + ...
			     	 ti * (zzd(1,1,1) + zzd(1,2,2) + zzd(2,1,1) + zzd(2,2,2)));

	options = sdpsettings('solver','mosek','verbose',0);
	sol = optimize(constraints,obj,options); % solve the SDP

	%% compute quasi-relative entropy
	f = wi/(ti*log(2)); % prefactor
	Hi = f * (1 + value(obj)); % quasi-relative entropy
end