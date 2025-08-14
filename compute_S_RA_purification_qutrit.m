function S_RA = compute_S_RA_purification_qutrit(p, R_pure, K, n)
    k = size(R_pure, 3);
    d_in = size(R_pure, 1);
    d_out = size(K{1}, 1);
    tol = 1e-10;

    psi_vectors = cell(k, 1);
    for i = 1:k
        [V, D] = eig(R_pure(:,:,i));
        [~, idx] = max(diag(D));
        psi_vectors{i} = V(:, idx);
    end

    sigma_ij = cell(k, k);
    for i = 1:k
        for j = 1:k
            op_ij = zeros(d_out, d_out);
            for l = 1:numel(K)
                op_ij = op_ij + K{l} * (psi_vectors{i} * psi_vectors{j}') * K{l}';
            end
            sigma_ij{i,j} = op_ij;
        end
    end

    partitions = generate_partitions(n, d_out);
    filename = sprintf('dimWlambda/dim_W_lambda_n%d_d%d.mat', n, d_out); 
    data = load(filename); 
    dim_S_lambda = double(data.dim_W_lambda);
    
    E_matrices_filename = sprintf('E_matrices/E_matrices_n%d_d%d.mat', n, d_out); 
    load(E_matrices_filename, 'E_matrices');
    
    num_partitions = size(partitions, 1);
    c_lambda_list = zeros(num_partitions, 1);
    entropy_sum = 0; 
    
    for part_idx = 1:num_partitions
        E = E_matrices{part_idx};
        q_lambda_dim = size(E.E11, 1);
        
        q_sigma_n = cell(k, k);
        parfor i = 1:k
            for j = 1:k
                A_ij = logm(sigma_ij{i,j});
                
                E13 = E.E12 * E.E23 - E.E23 * E.E12;
                E31 = E.E32 * E.E21 - E.E21 * E.E32;

                rep_ij = A_ij(1,1)*E.E11 + A_ij(1,2)*E.E12 + A_ij(1,3)*E13 + ...
                         A_ij(2,1)*E.E21 + A_ij(2,2)*E.E22 + A_ij(2,3)*E.E23 + ...
                         A_ij(3,1)*E31   + A_ij(3,2)*E.E32 + A_ij(3,3)*E.E33;
                
                q_sigma_n{i,j} = fastexpm(rep_ij);
            end
        end
        
        Q_lambda = zeros(k * q_lambda_dim);
        Q_lambda(1:q_lambda_dim, 1:q_lambda_dim) = p(1) * q_sigma_n{1,1};
        Q_lambda(1:q_lambda_dim, q_lambda_dim+1:end) = sqrt(p(1)*p(2)) * q_sigma_n{1,2};
        Q_lambda(q_lambda_dim+1:end, 1:q_lambda_dim) = sqrt(p(1)*p(2)) * q_sigma_n{2,1};
        Q_lambda(q_lambda_dim+1:end, q_lambda_dim+1:end) = p(2) * q_sigma_n{2,2};
        
        Q_lambda = (Q_lambda + Q_lambda') / 2;
        if norm(Q_lambda - Q_lambda', 'fro') > tol 
            warning('Q_lambda not Hermitian, diff=%.2e', norm(Q_lambda-Q_lambda','fro')); 
        end
        
        c_lambda_unweighted = trace(Q_lambda);
        
        if c_lambda_unweighted * dim_S_lambda(part_idx) > 1e-12
            c_lambda = c_lambda_unweighted * dim_S_lambda(part_idx);
            c_lambda_list(part_idx) = c_lambda;
            
            omega_lambda = Q_lambda / c_lambda_unweighted;
            S_omega_lambda = VNent(omega_lambda);
            log_dim_S = log2(dim_S_lambda(part_idx));
            
            entropy_sum = entropy_sum + c_lambda * (S_omega_lambda + log_dim_S);
        end
    end
    
    total_trace = sum(c_lambda_list);
    c_lambda_dist = c_lambda_list / total_trace;
    H_c = -sum(c_lambda_dist(c_lambda_dist > 0) .* log2(c_lambda_dist(c_lambda_dist > 0)));
    
    S_RA = H_c + entropy_sum;
end
