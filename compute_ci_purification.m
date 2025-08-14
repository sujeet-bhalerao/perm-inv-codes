function CI = compute_ci_purification(p, R_pure, K, n)
    d = size(R_pure, 1);
    k = size(R_pure, 3);


    rhoB_list = cell(1, k);
    for i = 1:k
        rhoB = zeros(d);
        for l = 1:numel(K)
            rhoB = rhoB + K{l} * R_pure(:,:,i) * K{l}';
        end
        rhoB_list{i} = rhoB;
    end
    S_B = real(compute_entropy_d2(n, d, k, rhoB_list, p));

    
    S_RA = real(compute_S_RA_GL2(p, R_pure, K, n));
    

    CI = real(S_B - S_RA) / n;
end
