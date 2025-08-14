function S_sigma_n = compute_entropy_d3(n, d, k, rho_list, p)
  

partitions = generate_partitions(n, d);
num_partitions = size(partitions, 1);

filename = sprintf('dimWlambda/dim_W_lambda_n%d_d%d.mat', n, d);
data = load(filename);
dim_W_lambda = double(data.dim_W_lambda);

E_matrices_filename = sprintf('E_matrices/E_matrices_n%d_d%d.mat', n, d);
if exist(E_matrices_filename, 'file')
    load(E_matrices_filename, 'E_matrices');
else

    E_matrices = generate_E_matrices_d3(n, d, partitions);
    save(E_matrices_filename, 'E_matrices');
   
end

results = cell(k, 1);

parfor rho_idx = 1:k
    rho = rho_list{rho_idx};
    A = q_logm(rho);

    result = struct();

    for part_idx = 1:num_partitions
        lambda = partitions(part_idx, :);
        lambda = lambda(lambda > 0);

        E = E_matrices{part_idx};

        rep = A(1,1) * E.E11 + A(1,2) * E.E12 + A(1,3) * (E.E12 * E.E23 - E.E23 * E.E12) + ...
              A(2,1) * E.E21 + A(2,2) * E.E22 + A(2,3) * E.E23 + ...
              A(3,2) * E.E32 + A(3,3) * E.E33 + A(3,1) * (E.E32 * E.E21 - E.E21 * E.E32);

        rep_exp = fastexpm(rep);

        result(part_idx).partition = lambda;
        result(part_idx).rep_exp = rep_exp;
    end

    results{rho_idx} = result;
end
num_partitions = size(partitions, 1);
c_lambda_values = zeros(num_partitions, 1);
entropy_contributions = zeros(num_partitions, 1);

parfor part_idx = 1:num_partitions
    rep_exp_size = size(results{1}(part_idx).rep_exp);
    q_lambda_bar_current = zeros(rep_exp_size);
    
    for rho_idx = 1:k
        p_i = p(rho_idx);
        rep_exp = results{rho_idx}(part_idx).rep_exp;
        q_lambda_bar_current = q_lambda_bar_current + p_i * rep_exp;
    end
   
    q_lambda_bar_current = q_lambda_bar_current * dim_W_lambda(part_idx);

    c_lambda = real(trace(q_lambda_bar_current));

    if c_lambda > 0
        c_lambda_values(part_idx) = c_lambda;
        q_normalized = q_lambda_bar_current / c_lambda;
 
        eigenvals = eig(q_normalized);
        eigenvals = real(eigenvals);
        positive_eigenvals = eigenvals(eigenvals > 0);
        
        if ~isempty(positive_eigenvals)
            normalized_eigenvals = positive_eigenvals / sum(positive_eigenvals);
            S_q_lambda = -sum(normalized_eigenvals .* log2(normalized_eigenvals));
            log_dim_W_lambda = log2(dim_W_lambda(part_idx));
            entropy_contributions(part_idx) = c_lambda * (S_q_lambda + log_dim_W_lambda);
        end
    end
end

c_lambda_list = c_lambda_values(c_lambda_values > 0);
entropy_sum = sum(entropy_contributions);

if isempty(c_lambda_list)
    S_sigma_n = entropy_sum;
else
    H_c = -sum(c_lambda_list .* log2(c_lambda_list));
    S_sigma_n = H_c + entropy_sum;
end

end



