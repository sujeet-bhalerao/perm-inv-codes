function S_sigma_n = compute_entropy_d2(n, d, k, rho_list, p)
 
    if length(rho_list) ~= k || length(p) ~= k
        error('The number of rhos and probabilities must match k');
    end
    
    if abs(sum(p) - 1) > 1e-10
        error('Probabilities must sum to 1');
    end
    
  
    partitions = generate_partitions(n, d);
    num_partitions = size(partitions, 1);
    
  
    filename = sprintf('dimWlambda/dim_W_lambda_n%d_d%d.mat', n, d);
    data = load(filename);
    dim_W_lambda = double(data.dim_W_lambda);
    
  
    E_matrices_filename = sprintf('E_matrices/E_matrices_n%d_d%d.mat', n, d);
    if ~exist(E_matrices_filename, 'file')
        error('E matrices file not found. Please generate it first.');
    end
    load(E_matrices_filename, 'E_matrices');

  
    results = cell(k, 1);
    for rho_idx = 1:k
        rho = rho_list{rho_idx};
        A = q_logm(rho);
        
        result = struct();
        for part_idx = 1:num_partitions
            lambda = partitions(part_idx, :);
            lambda = lambda(lambda > 0);
            
            E = E_matrices{part_idx};
            
            rep = A(1,1) * E.E11 + A(1,2) * E.E12 + A(2,1) * E.E21 + A(2,2) * E.E22;
            
            rep_exp = expm(rep);
            
            result(part_idx).partition = lambda;
            result(part_idx).rep_exp = rep_exp;
        end
        
        results{rho_idx} = result;
    end
    

    q_lambda_bar = cell(num_partitions, 1);
    c_lambda_list = [];
    entropy_sum = 0;
    
    for part_idx = 1:num_partitions
        rep_exp_size = size(results{1}(part_idx).rep_exp);
        q_lambda_bar{part_idx} = zeros(rep_exp_size);
        
        for rho_idx = 1:k
            p_i = p(rho_idx);
            rep_exp = results{rho_idx}(part_idx).rep_exp;
            q_lambda_bar{part_idx} = q_lambda_bar{part_idx} + p_i * rep_exp;
        end
        
        q_lambda_bar{part_idx} = q_lambda_bar{part_idx} * dim_W_lambda(part_idx);
        
        c_lambda = trace(q_lambda_bar{part_idx});
        
        if c_lambda > 0
            c_lambda_list(end+1) = c_lambda;
            
            q_normalized = q_lambda_bar{part_idx} / c_lambda;
            
            eigenvals = eig(q_normalized);
            eigenvals = real(eigenvals);
            eigenvals = eigenvals(eigenvals > 0);
            eigenvals = eigenvals / sum(eigenvals);
            
            S_q_lambda = -sum(eigenvals .* log2(eigenvals));
            log_dim_W_lambda = log2(dim_W_lambda(part_idx));
            
            entropy_sum = entropy_sum + c_lambda * (S_q_lambda + log_dim_W_lambda);
        end
    end
    
   
    c_lambda_list = c_lambda_list / sum(c_lambda_list);
    
    
    H_c = -sum(c_lambda_list .* log2(c_lambda_list));
    
  
    S_sigma_n = H_c + entropy_sum;
end

function partitions = generate_partitions(n, d)
    partitions = [];
    for i = n:-1:0
        j = n - i;
        if j <= i && length([i, j]) <= d
            partitions = [partitions; [i, j]];
        end
    end
end



