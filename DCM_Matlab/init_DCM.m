function[found, it, time] = init_DCM(M1, M2, k, del, start)
	% Takes pre-prepared Matrix (standardized and optionally Quantile Normalized). Looks for high correlation in M1 and low/negative correlation in M2, finds group of size k.
    
   if exist('del', 'var')
		% Record dimension
		real_p = size(M1, 1);
		idcs = (1:real_p);
        idcs(del) = [];

		% Remove ignored rows
		M1(del, :) = [];
		M2(del, :) = [];
   else
		% Record dimension
		p = size(M1, 1);
		idcs = 1:p;
   end

	% Start with prespecified row set or with a random k rows.
	if exist('start', 'var')
        print('hi')
		A = start;
        k = length(A);
    else
		A = randsample(length(idcs), k);
    end

    % Sizes
    [p, n1] = size(M1);
    n2 = size(M2, 2);
	
	% Save starting point
	orig_A = A;

	% Make list of all indices, and of those not in seed set
	P = 1:p;
	notA = P(~ismember(P, A));

	% Initialize for loop
	done = false;
	it = 0;

	% Find correlations of all genes with only genes in A
	cross_1 = round(M1*M1(A,:).', 10);
	cross_1(cross_1 == 1) = 0; % So that transformed value will be 0 and variance won't contribute to sum
	cross_1 = fisher(cross_1)*sqrt(n1 - 3); % Fisher transform
	cross_2 = round(M2*M2(A,:).', 10);
	cross_2(cross_2 == 1) = 0;
	cross_2 = fisher(cross_2)*sqrt(n2 - 3);

	% Note: resulting matrices have rows in order of data indices, columns correspond to values of A in the order listed by A

	% Rowsums represent total of pairwise correlations between [row] and A	
	rows_1 = sum(cross_1, 2);
	rows_2 = sum(cross_2, 2);

	d = (p-k)*k; % Number of in/out swap combos
	
	% Count iterations
	it = 0;
	
	% Time it
	tic;
	
	% Iterate until convergence
	while ~done
    
      	% Rowsum Differences
		diffs12 = rows_1 - rows_2;
			
		% Gain due to o is getting back the contribution of o from h, loss due to o is corrs for s
		% Similarly for in
        for i = 1:k
			temp = diffs12(notA) + cross_2(notA, i) - cross_1(notA, i);
            [a,b] = max(temp);
            bestAs(i, 1:2) = [b,a];
        end

		bestAs(:,2) = bestAs(:,2) - diffs12(A);
		[m, best_out] = max(bestAs(:,2));
		best_in = bestAs(best_out,1);

		if bestAs(best_out,2) <= 0
			done = true;
        else
			% Find in and out data labels
			
			out = A(best_out);
			inn = notA(best_in);
		
			% Switch 'out' and 'in' indices from their lists
			A(best_out) = inn;
			notA(best_in) = out;
		
			% Find correlation vector for 'in'
			new_1 = round(M1*M1(inn,:).', 10);
			new_1(inn) = 0; % So Fisher isn't inf
			new_1 = fisher(new_1)*sqrt(n1 - 3);
			% Edit rowsums to reflect inclusion of new index, exclusion of old
			rows_1 = rows_1 - cross_1(:, best_out) + new_1;
			% Replace relevant column of corr matrix with new correlations
			cross_1(:, best_out) = new_1;
		
			% Same thing for Group 2
			new_2 = round(M2*M2(inn, :).', 10);
			new_2(inn) = 0;
			new_2 = fisher(new_2)*sqrt(n2 - 3);
			rows_2 = rows_2 - cross_2(:, best_out) + new_2;
			cross_2(:, best_out) = new_2;
			
			% Increase iteration count
			it = it + 1;
        end
    end
	time = toc;
	
	% Translate back to real indices
	found = idcs(A).';
	
end
