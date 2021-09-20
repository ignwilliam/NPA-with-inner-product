function [ops] = GeneratePartyOps(party,n_in,n_out)

	in_cell = cell(n_in,n_out);
	out_cell = cell(n_in,n_out);
	for i = 1:n_in
		for j = 1:n_out
			in_cell{i,j} = i-1;
			out_cell{i,j} = j-1;
		end
	end
	
	switch party
		case 'A'
			ops = struct('status','1','as',in_cell,'ao',out_cell,'bs',[],'bo',[],...
					     'cs',[],'co',[],'cdagger',[]);
		case 'B'
			ops = struct('status','1','as',[],'ao',[],'bs',in_cell,'bo',out_cell,...
						 'cs',[],'co',[],'cdagger',[]);
		otherwise
			ops = struct('status','1','as',[],'ao',[],'bs',[],'bo',[],...
						 'cs',in_cell,'co',out_cell,'cdagger',0);
	end
end