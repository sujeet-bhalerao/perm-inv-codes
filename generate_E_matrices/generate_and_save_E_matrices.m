function generate_and_save_E_matrices(n, d)
    partitions = generate_partitions(n, d);
    E_matrices = generate_E_matrices_sparse(n, d, partitions);

    if ~exist('E_matrices_sparse', 'dir')
        mkdir('E_matrices_sparse');
    end
    
    filename = sprintf('E_matrices_sparse/E_matrices_sparse_n%d_d%d.mat', n, d);
    save(filename, 'E_matrices', '-v7.3');  
    fprintf('Sparse E matrices saved to %s\n', filename);
end
