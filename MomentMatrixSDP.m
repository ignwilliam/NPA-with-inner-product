function constraints = MomentMatrixSDP(v, Gamma, constraints, lambda, idx_matrix)

	L = size(idx_matrix,1);
	num_states = size(lambda,1);
	d = L * num_states;

	for i = 1:num_states
		for j = 1:num_states
			for k = 1:L
				for l = 1:L
					row = (i-1)*L + k;
					col = (j-1)*L + l;
					
					switch idx_matrix(k,l)
						case 0
							constraints = [constraints, Gamma(row,col) == 0];
						case 1
							constraints = [constraints, Gamma(row,col) == lambda(i,j)];
						otherwise
							constraints = [constraints, Gamma(row,col) == v(idx_matrix(k,l),i,j)];
					end
				end
			end
		end
	end

	Gamma = (Gamma + Gamma')/2;
	constraints = [constraints, Gamma >= 0];
end