function E_matrices = generate_E_matrices_d3(n, d, partitions)
    if d ~= 3
        error('This function is for d=3');
    end

    num_partitions = size(partitions, 1);
    E_matrices = cell(num_partitions, 1);
    
    for part_idx = 1:num_partitions
        lambda = partitions(part_idx, :);
        lambda_active = lambda(lambda > 0);
        
        E = struct();
        
        E.E11 = sparse(E_kk_rep(n, d, lambda_active, 1));
        E.E22 = sparse(E_kk_rep(n, d, lambda_active, 2));
        E.E33 = sparse(E_kk_rep(n, d, lambda_active, 3));

        E.E12 = sparse(E_k_kp1_rep2_normalized(n, d, lambda_active, 1));
        E.E21 = sparse(E_kp1_k_rep2_normalized(n, d, lambda_active, 1));

        E.E23 = sparse(E_k_kp1_rep2_normalized(n, d, lambda_active, 2));
        E.E32 = sparse(E_kp1_k_rep2_normalized(n, d, lambda_active, 2));
        
        E_matrices{part_idx} = E;
    end
end
