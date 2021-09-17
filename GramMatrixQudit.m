function lambda = GramMatrixQudit(states)

	num_states = size(states,1);
	lambda = zeros(num_states,num_states);
	
	for i = 1:num_states
		for j = 1:num_states
			lambda(i,j) = states(i,:) * states(j,:)';
		end
	end

end