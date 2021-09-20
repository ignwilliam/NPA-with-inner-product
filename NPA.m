classdef NPA
	properties
		gram_input
		level
		io_config
		monomials
		isComplex

		moment_matrix
		ref_table
		index_matrix

		variables
		Gamma
		npa_constraints
	end

	methods
		%% constructor for the NPA with inner-product class
		function sdp = NPA(gram_input, io_config, level, extra_monomials, isComplex)
			% Gram matrix of the input states, input/output configuration, and relaxation level
			sdp.gram_input = gram_input;
			sdp.io_config = io_config;
			sdp.level = level;
			sdp.isComplex = isComplex;

			% generate monomials
			sdp.monomials = GenerateOps(io_config, level);
			sdp.monomials = [sdp.monomials, extra_monomials];

			% generate the moment matrix
			[sdp.moment_matrix, sdp.ref_table, sdp.index_matrix] = GenerateMomentMatrix(sdp.monomials);

			% calculate sizes
			num_mono = length(sdp.monomials);
			num_states = size(gram_input,1);
			num_var = length(sdp.ref_table);

			% generate SDP variables
			if isComplex
				sdp.variables = sdpvar(num_var,num_states,num_states, 'full', 'complex');
				sdp.Gamma = sdpvar(num_mono*num_states, num_mono*num_states, 'full', 'complex');
			else
				sdp.variables = sdpvar(num_var,num_states,num_states, 'full');
				sdp.Gamma = sdpvar(num_mono*num_states, num_mono*num_states, 'full');
			end
			
			sdp.npa_constraints = MomentMatrixSDP(sdp.variables, sdp.Gamma, [], sdp.gram_input, sdp.index_matrix);
		end

		function vars = Operator2Variable(sdp, ops, bra, ket)
			
			% calculate dimension of the matrix of operators, ops
			dims = size(ops);
			
			% flatten ops
			ops = reshape(ops,[1,numel(ops)]);

			% convert to variables
			for i = 1:length(ops)
				vars(i) = sdp.variables(IndexOps(ops(i),sdp.ref_table),bra,ket);
			end

			% reshaping the output
			vars = reshape(vars,dims);
		end



	end

end