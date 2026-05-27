function Ekp1k = E_kp1_k_rep2_normalized(n, d, lambda, k)
    % Function to compute the normalized matrix representation of E_{k+1,k}
    % Inputs:
    %   n      - Total sum of parts
    %   d      - Maximum number of parts allowed 
    %   lambda - Partition of n into at most d parts (row vector)
    %   k      - The specific k value (integer from 1 to d-1)
    % Output:
    %   Ekp1k  - The normalized matrix representation of E_{k+1,k}
    
    lambda = lambda(:)'; 
    if length(lambda) > d
        error('The partition lambda has more than d parts.');
    end
    lambda_padded = [lambda, zeros(1, d - length(lambda))];
    

    if k < 1 || k >= d || floor(k) ~= k
        error('k must be an integer between 1 and d-1.');
    end
    
    filename = sprintf('gt_patterns/gt_patterns_n%d_d%d.mat', n, d);
    if exist(filename, 'file') ~= 2
        error('File %s does not exist. Generate the GT patterns and save them as a .mat file.', filename);
    end
    data = load(filename);
    gt_patterns = data.gt_patterns;
    

    [num_patterns, ~, ~] = size(gt_patterns);
    

    matching_patterns = {};
    for i = 1:num_patterns
        pattern = squeeze(gt_patterns(i, :, :)); 
        top_row = pattern(1, :);
        if isequal(top_row, lambda_padded)
            matching_patterns{end+1} = pattern;
        end
    end
    
    num_matching_patterns = length(matching_patterns);
    if num_matching_patterns == 0
        error('No GT patterns found with the given partition lambda.');
    end
    
    pattern_to_index = containers.Map();
    for idx = 1:num_matching_patterns
        pattern = matching_patterns{idx};
        key = pattern_to_key(pattern);
        pattern_to_index(key) = idx;
    end
    
    pattern_norms = zeros(num_matching_patterns, 1);
    for idx = 1:num_matching_patterns
        pattern = matching_patterns{idx};
        pattern_norms(idx) = compute_inner_product_alld(pattern);
    end
    
    row_indices = [];
    col_indices = [];
    values = [];
    

    for idx = 1:num_matching_patterns
        pattern = matching_patterns{idx};
        norm_current = pattern_norms(idx);
        

        for i = 1:k
            lambda_ki = pattern(d - k + 1, i);
            l_ki = double(lambda_ki - i + 1); 
            
            numerator = 1.0;  
            for j = 1:(k - 1)

                row_km1_idx = d - (k - 1) + 1;
                if row_km1_idx >= 1 && row_km1_idx <= d
                    if j <= d - row_km1_idx + 1
                        lambda_km1j = pattern(row_km1_idx, j);
                        l_km1j = double(lambda_km1j - j + 1);
                    else
                        l_km1j = 0.0;
                    end
                else
                    l_km1j = 0.0;
                end
                numerator = numerator * (l_ki - l_km1j);
            end
            

            denominator = 1.0;
            for m = 1:k
                if m ~= i
                    lambda_km = pattern(d - k + 1, m);
                    l_km = double(lambda_km - m + 1);
                    denominator = denominator * (l_ki - l_km);
                end
            end
            
            if denominator == 0
                coeff = 0.0;
            else
                coeff = numerator / denominator;
            end
            
            if coeff ~= 0.0
      
                new_pattern = pattern;
                new_pattern(d - k + 1, i) = new_pattern(d - k + 1, i) - 1;
                
           
                if is_valid_gt_pattern(new_pattern)
          
                    key_new = pattern_to_key(new_pattern);
                    if isKey(pattern_to_index, key_new)
                        idx_new = pattern_to_index(key_new);
                        norm_new = pattern_norms(idx_new);
                        
            
                        adjustment_factor = sqrt(norm_new / norm_current);
                        

                        coeff_adjusted = coeff * adjustment_factor;
                        

                        row_indices(end+1, 1) = idx_new;
                        col_indices(end+1, 1) = idx;
                        values(end+1, 1) = coeff_adjusted;
                    end

                end
            end
        end
    end
    
    Ekp1k = sparse(row_indices, col_indices, values, num_matching_patterns, num_matching_patterns);
end
