function Ekk = E_kk_rep(n, d, lambda, k)
    % Function to compute the E_kk matrix for a given partition lambda and specific k
    % Inputs:
    %   n      - Total sum of parts
    %   d      - Maximum number of parts allowed
    %   lambda - Partition of n into at most d parts (row vector)
    %   k      - The specific k value (integer from 1 to d)
    % Output:
    %   Ekk    - The E_kk matrix (diagonal matrix)


    lambda = lambda(:)';
    if length(lambda) > d
        error('The partition lambda has more than d parts.');
    end
    lambda_padded = [lambda, zeros(1, d - length(lambda))];

    if k < 1 || k > d || floor(k) ~= k
        error('k must be an integer between 1 and d.');
    end


    filename = sprintf('gt_patterns/gt_patterns_n%d_d%d.mat', n, d);
    

    if exist(filename, 'file') ~= 2
        error('File %s does not exist. Generate the GT patterns and save them as a .mat file.', filename);
    end
    
    % disp(['Loading GT patterns from file: ', filename]);
    data = load(filename);
    gt_patterns = data.gt_patterns;
    

    [num_patterns, num_rows, num_columns] = size(gt_patterns);
    % disp(['Number of patterns: ', num2str(num_patterns)]);
    % disp(['Rows: ', num2str(num_rows), ', Columns: ', num2str(num_columns)]);
    
    assert(num_rows == d, 'Row size mismatch! Expected: %d, but got: %d', d, num_rows);
    assert(num_columns == d, 'Column size mismatch! Expected: %d, but got: %d', d, num_columns);



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

    diag_entries = zeros(num_matching_patterns, 1);

    
    for i = 1:num_matching_patterns
        pattern = matching_patterns{i}; 

        % Indices for the rows
        row1_idx = d - k + 1;
        row2_idx = d - k + 2;

        % Sum of entries in (d - k + 1)th row
        if row1_idx >= 1 && row1_idx <= d
            sum_row1 = sum(pattern(row1_idx, :));
        else
            sum_row1 = 0;
        end

        % Sum of entries in (d - k + 2)th row
        if row2_idx >= 1 && row2_idx <= d
            sum_row2 = sum(pattern(row2_idx, :));
        else
            sum_row2 = 0;
        end

        diag_entries(i) = sum_row1 - sum_row2;
    end

    Ekk = diag(diag_entries);
end
