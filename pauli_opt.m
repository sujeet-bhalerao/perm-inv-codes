
d = 2; n = 9; 
k = 2;
nx = 4*k ;
p = 0.2271;
px = p/2;
py = 0;
pz = p/2;



timestamp = datestr(now, 'yyyymmdd_HHMMSS');
output_dir = fullfile(pwd, ['optimization_results/pauli_', timestamp]);
mkdir(output_dir);


param_file = fullfile(output_dir, 'parameters.txt');
fid = fopen(param_file, 'w');
fprintf(fid, 'Parameters:\n');
fprintf(fid, 'd = %d\n', d);
fprintf(fid, 'n = %d\n', n);
fprintf(fid, 'k = %d\n', k);
fprintf(fid, 'p = %.4f\n', p);
fprintf(fid, 'px = %.4f\n', px);
fprintf(fid, 'py = %.4f\n', py);
fprintf(fid, 'pz = %.4f\n', pz);
fclose(fid);

% Kraus operators for 2-Pauli channel
K0 = sqrt(1 - px - py - pz) * [1, 0; 0, 1];
K1 = sqrt(px) * [0, 1; 1, 0];
K2 = sqrt(py) * [0, -i; i, 0];
K3 = sqrt(pz) * [1, 0; 0, -1];
K = {K0, K1, K3};


Kc = comp_kraus(K,[2,2]);

obj = @(x) -compute_ci_symmetries_opt(x,d,k,K,Kc,n);

maxIterations = 250;
swarmSize = 200;

opt = optimoptions('particleswarm', ...
    'Display', 'iter', ...
    'SwarmSize', swarmSize, ...
    'MaxIterations', maxIterations, ...
    'MaxStallIterations', maxIterations, ...
    'FunctionTolerance',1e-10,...
    'OutputFcn', @(optimValues,state) save_iteration(optimValues,state,output_dir), ...
    'UseParallel', true);  


[x,f] = particleswarm(obj,nx,[],[],opt);


[p,R] = get_states_bloch(x,d,k);


save(fullfile(output_dir, 'final_results.mat'), 'x', 'f', 'p', 'R');


result_file = fullfile(output_dir, 'final_results.txt');
fid = fopen(result_file, 'w');
fprintf(fid, 'Final Results:\n');
fprintf(fid, 'Objective value: %.6f\n\n', f);
fprintf(fid, 'Optimized probabilities:\n');
fprintf(fid, '%f\n', p);
fprintf(fid, '\nOptimized states:\n');
for i = 1:k
    fprintf(fid, 'R%d:\n', i);
    fprintf(fid, '%f %f\n%f %f\n\n', R(:,:,i)');
end
fclose(fid);


function stop = save_iteration(optimValues,state,output_dir)
    stop = false;
    iteration_file = fullfile(output_dir, 'iteration_data.csv');
    
    if strcmp(state, 'init')
       
        fid = fopen(iteration_file, 'w');
        fprintf(fid, 'Iteration,Best Objective,Mean Objective\n');
        fclose(fid);
    end
    
    if strcmp(state, 'iter')
        
        fid = fopen(iteration_file, 'a');
        fprintf(fid, '%d,%.6f,%.6f\n', optimValues.iteration, optimValues.bestfval, mean(optimValues.swarmfvals));
        fclose(fid);
    end
end
