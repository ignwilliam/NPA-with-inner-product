function h = h2(p) %binary entropy

h = zeros(1,length(p));

% check if p = 0 or p = 1
i = (p ~= 0) & (p ~= 1);

% binary entropy
h(i) = - p(i) .* log2(p(i)) - (1 - p(i)) .* log2(1 - p(i));