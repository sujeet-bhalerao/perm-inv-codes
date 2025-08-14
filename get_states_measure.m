function [p,R] = get_states_measure(x,d,k)
% x ... variable vector
% d ... channel input dimension
% k ... number of states in convex combination
R = zeros(d,d,k);
p = zeros(1,k);
N = k*d^2;
psi(:,1) = x(1:N)+1i*x(N+1:end);
psi = psi/norm(psi);
Ik = eye(k);
for j=1:k
    phi = kron(Ik(:,j)',eye(d^2))*psi;
    p(j) = phi'*phi;
    phi = phi/norm(phi);
    R(:,:,j) = TrX(phi,2,[d,d]);
end
end