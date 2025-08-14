function res = compute_ci_symmetries_opt(x, d, k, K, Kc, n, precomp)
    [p, R] = get_states_bloch(x, d, k);
    epsilon = 1e-12;
    
    
    I_d = eye(d);


    mB = length(K);
    dB = size(K{1}, 1);
    mE = length(Kc);
    dE = size(Kc{1}, 1);
    
    rho_list_B = cell(1, k);
    rho_list_E = cell(1, k);
    
    for j = 1:k        
        R(:,:,j) = (1 - epsilon) * R(:,:,j) + epsilon * I_d / d;
        opB = zeros(dB, dB);
        for l = 1:mB
            term = K{l} * R(:,:,j) * K{l}';
            opB = opB + term;
        end
        rho_list_B{j} = opB;
        
        opE = zeros(dE, dE);
        for l = 1:mE
            term = Kc{l} * R(:,:,j) * Kc{l}';
            opE = opE + term;
        end
        rho_list_E{j} = opE;
    end

    dimB = size(rho_list_B{1}, 1);
    if dimB == 2
        S_sigmaB = compute_entropy_d2(n, dimB, k, rho_list_B, p);
    elseif dimB == 3
        S_sigmaB = compute_entropy_d3(n, dimB, k, rho_list_B, p);
    elseif dimB == 4
        S_sigmaB = compute_entropy_d4(n, dimB, k, rho_list_B, p, precomp);
    else
        error('No entropy routine available for channel B dimension %d', dimB);
    end

    dimE = size(rho_list_E{1}, 1);
    if dimE == 2
        S_sigmaE = compute_entropy_d2(n, dimE, k, rho_list_E, p);
    elseif dimE == 3
        S_sigmaE = compute_entropy_d3(n, dimE, k, rho_list_E, p);
    elseif dimE == 4
        S_sigmaE = compute_entropy_d4(n, dimE, k, rho_list_E, p, precomp);
    else
        error('No entropy routine available for channel E dimension %d', dimE);
    end

    res = real(S_sigmaB - S_sigmaE) / n;
end
