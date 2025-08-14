
function S_RA = compute_S_RA_log(p, R_pure, K, n)
    k=size(R_pure,3);
    d=size(R_pure,1);
    tol=1e-10;

    psi_vectors=cell(k,1);
    for i=1:k,[V,~]=eig(R_pure(:,:,i));psi_vectors{i}=V(:,end);end
    sigma_ij=cell(k,k);
    for i=1:k
        for j=1:k
            op_ij=zeros(d);
            for l=1:numel(K)
                op_ij=op_ij+K{l}*(psi_vectors{i}*psi_vectors{j}')*K{l}';
            end
            sigma_ij{i,j}=op_ij;
        end
    end
    partitions=generate_partitions(n,d);
    filename=sprintf('dimWlambda/dim_W_lambda_n%d_d%d.mat',n,d);
    data=load(filename);
    dim_S_lambda=double(data.dim_W_lambda);
    E_matrices_filename=sprintf('E_matrices/E_matrices_n%d_d%d.mat',n,d);
    load(E_matrices_filename,'E_matrices');

    num_partitions=size(partitions,1);
    c_lambda_list=zeros(num_partitions,1);
    entropy_sum = 0; 
    
    parfor part_idx=1:num_partitions
        E=E_matrices{part_idx};
        q_lambda_dim=size(E.E11,1);
        q_sigma_n=cell(k,k);
        for i=1:k
            for j=1:k
                A_ij=logm(sigma_ij{i,j});
                rep_ij=A_ij(1,1)*E.E11+A_ij(1,2)*E.E12+A_ij(2,1)*E.E21+A_ij(2,2)*E.E22;
                q_sigma_n{i,j}=expm(rep_ij);
            end
        end
        Q_lambda=zeros(k*q_lambda_dim);
        Q_lambda(1:q_lambda_dim,1:q_lambda_dim)=p(1)*q_sigma_n{1,1};
        Q_lambda(1:q_lambda_dim,q_lambda_dim+1:end)=sqrt(p(1)*p(2))*q_sigma_n{1,2};
        Q_lambda(q_lambda_dim+1:end,1:q_lambda_dim)=sqrt(p(1)*p(2))*q_sigma_n{2,1};
        Q_lambda(q_lambda_dim+1:end,q_lambda_dim+1:end)=p(2)*q_sigma_n{2,2};
        Q_lambda = (Q_lambda + Q_lambda') / 2;

        herm_diff=norm(Q_lambda-Q_lambda','fro')
        if herm_diff > tol 
        
            warning('Q_lambda not Hermitian, diff=%.2e',herm_diff); 
        end
        
        c_lambda_unweighted = trace(Q_lambda);
        c_lambda = c_lambda_unweighted * dim_S_lambda(part_idx);
        
        if c_lambda > 0
            c_lambda_list(part_idx) = c_lambda;
            
            omega_lambda = Q_lambda / c_lambda_unweighted;
            S_omega_lambda = VNent(omega_lambda);
            log_dim_S = log2(dim_S_lambda(part_idx));
            
   
            entropy_sum = entropy_sum + c_lambda * (S_omega_lambda + log_dim_S);
        end
    end
    
    total_trace = sum(c_lambda_list);
    c_lambda_dist = c_lambda_list / total_trace;
    H_c = -sum(c_lambda_dist(c_lambda_dist > 0) .* log2(c_lambda_dist(c_lambda_dist > 0)));
    
    S_RA = H_c + entropy_sum;
end
