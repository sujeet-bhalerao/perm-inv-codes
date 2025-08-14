
function logRho = q_logm(rho)
    [U, D] = eig(rho);
    D = diag(D);
    logD = zeros(size(D));
    logD(D > 0) = log(D(D > 0));
    logRho = U * diag(logD) * U';
end
