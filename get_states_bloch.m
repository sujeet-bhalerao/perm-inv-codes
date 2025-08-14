function [p,R] = get_states_bloch(x, d, k)
if d ~= 2
    error('get_states_bloch supports only qubit (d=2) states.');
end
% Total parameters: 4k = k probability params + 3k state parameters.
probParams = x(1:k);
stateParams = x(k+1:end);

% enforce non-negativity
p = abs(probParams(:));
% enforce normalization
p = p/sum(p);

R = zeros(2,2,k);
I = eye(2);
sigma_x = [0,1;1,0];
sigma_y = [0,-1i;1i,0];
sigma_z = [1,0;0,-1];
for j = 1:k
    idx = (j-1)*3;
    xx = stateParams(idx+1:idx+3);
    norm_val = norm(xx);
    r = xx/(1+norm_val);
    R(:,:,j) = (I + r(1)*sigma_x + r(2)*sigma_y + r(3)*sigma_z) / 2;
end

tol = 1e-10;
if abs(sum(p)-1) > tol || any(p < -tol)
    error('Probability vector is not valid. Sum = %g, negative elements exist.', sum(p));
end

for j = 1:k
    rj = R(:,:,j);
    if abs(trace(rj)-1) > tol
        error('Density matrix R(:,:, %d) does not have trace 1 (trace = %g).', j, trace(rj));
    end
    if norm(rj - rj','fro') > tol
        error('Density matrix R(:,:, %d) is not Hermitian.', j);
    end
    eigsR = eig(rj);
    if any(eigsR < -tol)
        error('Density matrix R(:,:, %d) is not PSD.', j);
    end
end
end
