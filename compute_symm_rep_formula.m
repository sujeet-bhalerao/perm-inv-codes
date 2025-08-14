function Sm = compute_symm_rep_formula(A, m)
    m = round(m);
    rep_dim = m + 1;
    Sm = zeros(rep_dim, rep_dim, 'like', A);
    a = A(1,1); b = A(1,2);
    c = A(2,1); d = A(2,2);

    for j_idx = 1:rep_dim
        for k_idx = 1:rep_dim
            j = j_idx - 1;
            k = k_idx - 1;
            
            p_min = max(0, k + j - m);
            p_max = min(k, j);
            
            summand = 0;
            for p = p_min:p_max
                term = nchoosek(j, p) * nchoosek(m - j, k - p) * ...
                       a^p * c^(j - p) * b^(k - p) * d^(m - j - (k - p));
                summand = summand + term;
            end
            
            norm_factor = sqrt(factorial(k)*factorial(m-k)) / ...
                          sqrt(factorial(j)*factorial(m-j));
                          
            Sm(k_idx, j_idx) = norm_factor * summand;
        end
    end
end
