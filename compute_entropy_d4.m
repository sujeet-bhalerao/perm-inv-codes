function S_sigma_n = compute_entropy_d4(n, d, k, rho_list, p, precomp)
    if length(rho_list) ~= k || length(p) ~= k
        error('The number of rhos and probabilities must match k');
    end
    if abs(sum(p) - 1) > 1e-10
        error('Probabilities must sum to 1');
    end

    partitions = precomp.partitions;
    num_partitions = size(partitions, 1);
    
    dim_W_lambda = precomp.dim_W_lambda;
    E_matrices = precomp.E_matrices;
    
    E_matrices_const = parallel.pool.Constant(E_matrices);
    partitions_const = parallel.pool.Constant(partitions);

    results = cell(k, 1);
    parfor rho_idx = 1:k  
        rho = rho_list{rho_idx};
        A = q_logm(rho);
        temp_result = repmat(struct('partition', [], 'rep_exp', []), num_partitions, 1);
        for part_idx = 1:num_partitions
            lambda = partitions_const.Value(part_idx, :); 
            lambda = lambda(lambda > 0);
            E = E_matrices_const.Value{part_idx}; 
            
            comm12_23 = E.E12 * E.E23 - E.E23 * E.E12;
            comm23_34 = E.E23 * E.E34 - E.E34 * E.E23;
            comm32_21 = E.E32 * E.E21 - E.E21 * E.E32;
            term14    = E.E12 * comm23_34 - comm23_34 * E.E12;
            term41    = E.E43 * comm32_21 - comm32_21 * E.E43;
            term42    = E.E43 * E.E32 - E.E32 * E.E43;
            
            rep = A(1,1) * E.E11 + A(2,2) * E.E22 + A(3,3) * E.E33 + A(4,4) * E.E44 + ...
                  A(1,2) * E.E12 + A(2,1) * E.E21 + ...
                  A(2,3) * E.E23 + A(3,2) * E.E32 + ...
                  A(3,4) * E.E34 + A(4,3) * E.E43 + ...
                  A(1,3) * comm12_23 + A(2,4) * comm23_34 + A(3,1) * comm32_21 + ...
                  A(1,4) * term14 + A(4,1) * term41 + A(4,2) * term42;
              
            rep_exp = fastexpm(rep);
            temp_result(part_idx).partition = lambda;
            temp_result(part_idx).rep_exp = rep_exp;
        end
        results{rho_idx} = temp_result;
    end

    local_c = zeros(num_partitions,1);
    local_entropy = zeros(num_partitions,1);
    
    parfor part_idx = 1:num_partitions
        rep_exp_size = size(results{1}(part_idx).rep_exp);
        q_local = zeros(rep_exp_size);
        for rho_idx = 1:k
            q_local = q_local + p(rho_idx) * results{rho_idx}(part_idx).rep_exp;
        end
        q_local = q_local * dim_W_lambda(part_idx);
        c_local = real(trace(q_local));
        if c_local > 0
            local_c(part_idx) = c_local;
            q_normalized = q_local / c_local;
            eigenvals = eig(q_normalized);
            eigenvals = real(eigenvals(real(eigenvals) > 0));
            eigenvals = eigenvals / sum(eigenvals);
            S_q_lambda = -sum(eigenvals .* log2(eigenvals));
            local_entropy(part_idx) = c_local * (S_q_lambda + log2(dim_W_lambda(part_idx)));
        else
            local_c(part_idx) = 0;
            local_entropy(part_idx) = 0;
        end
    end

    nonzero_c = local_c(local_c > 0);
    if isempty(nonzero_c)
        H_c = 0;
    else
        norm_c = nonzero_c / sum(nonzero_c);
        H_c = -sum(norm_c .* log2(norm_c));
    end
    entropy_sum = sum(local_entropy);
    S_sigma_n = H_c + entropy_sum;
end