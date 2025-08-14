function CI = compute_ci_symmetries(p,R_input,K,Kc,n,precomp) 
    if nargin<6, precomp = []; end
   
    k  = size(R_input,3);   
    mB = numel(K);
    mE = numel(Kc);


    if mB > 0
        d_out_B = size(K{1}, 1); 
    else
        d_out_B = size(R_input,1); 
    end
    
    rhoB = cell(1,k);
    rhoE = cell(1,k);

    for j=1:k
        opB = zeros(d_out_B); 
        for l=1:mB
            opB = opB + K{l}*R_input(:,:,j)*K{l}'; 
        end
        rhoB{j} = opB; 

        if mE > 0
            d_out_E = size(Kc{1}, 1);
        else
            d_out_E = size(R_input,1); 
             if isempty(Kc)
                warning('Kraus operator set Kc is empty.');
            end
        end
        opE = zeros(d_out_E);
        for l=1:mE
            opE = opE + Kc{l}*R_input(:,:,j)*Kc{l}';
        end
        rhoE{j} = opE; 
    end

    dimB = size(rhoB{1},1); 
    switch dimB
      case 2, S_B = compute_entropy_d2(n,dimB,k,rhoB,p);
      case 3, S_B = compute_entropy_d3(n,dimB,k,rhoB,p); 
      case 4, S_B = compute_entropy_d4(n,dimB,k,rhoB,p,precomp);
      otherwise, error('no routine for B dim %d',dimB)
    end

    dimE = size(rhoE{1},1);
    switch dimE
      case 2, S_E = compute_entropy_d2(n,dimE,k,rhoE,p);
      case 3, S_E = compute_entropy_d3(n,dimE,k,rhoE,p); 
      case 4, S_E = compute_entropy_d4(n,dimE,k,rhoE,p,precomp);
      otherwise, error('no routine for E dim %d',dimE)
    end

    CI = real(S_B - S_E)/n;
end
