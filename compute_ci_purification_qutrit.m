function CI = compute_ci_purification_qutrit(p, R_pure, K, n)
    d_out = size(K{1}, 1);
    k = size(R_pure, 3);

    rhoB_list = cell(1, k);
    for i = 1:k
        rhoB = zeros(d_out);
        for l = 1:numel(K)
            rhoB = rhoB + K{l} * R_pure(:,:,i) * K{l}';
        end
        rhoB_list{i} = rhoB;
    end
    
    S_B = compute_entropy_d3(n, d_out, k, rhoB_list, p);
    S_RA = compute_S_RA_purification_qutrit(p, R_pure, K, n);
    
    CI = real(S_B - S_RA) / n;
end
