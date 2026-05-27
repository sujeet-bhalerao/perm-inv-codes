function E_matrices = generate_E_matrices_sparse(n, d, partitions)
    num_partitions = size(partitions, 1);
    E_matrices = cell(num_partitions, 1);
    
    for part_idx = 1:num_partitions
        lambda = partitions(part_idx, :);
        lambda = lambda(lambda > 0);
        
        E = struct();
        E.E11 = sparse(E_kk_rep(n, d, lambda, 1));
        E.E22 = sparse(E_kk_rep(n, d, lambda, 2));
        E.E33 = sparse(E_kk_rep(n, d, lambda, 3));
        E.E44 = sparse(E_kk_rep(n, d, lambda, 4));
        
        E.E12 = sparse(E_k_kp1_rep2_normalized(n, d, lambda, 1));
        E.E21 = sparse(E_kp1_k_rep2_normalized(n, d, lambda, 1));
        
        E.E23 = sparse(E_k_kp1_rep2_normalized(n, d, lambda, 2));
        E.E32 = sparse(E_kp1_k_rep2_normalized(n, d, lambda, 2));
        
        E.E34 = sparse(E_k_kp1_rep2_normalized(n, d, lambda, 3));
        E.E43 = sparse(E_kp1_k_rep2_normalized(n, d, lambda, 3));
        
        E_matrices{part_idx} = E;
    end
end
