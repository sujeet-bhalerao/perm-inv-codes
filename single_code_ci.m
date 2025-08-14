
d = 2; n = 9; k = 2; 
p = [0.5, 0.5];
p_channel=0.2271; 
px=p_channel/2; 
pz=px; 
py=0;
K0=sqrt(1-px-py-pz)*eye(2); K1=sqrt(px)*[0 1; 1 0];K2 = sqrt(py) * [0, -i; i, 0]; K3=sqrt(pz)*[1 0; 0 -1];
K={K0, K1, K2, K3};

 

% code for 2-pauli channel
RR(:,:,1)=[0.179181326231154+0i,-0.029852749784765-0.382339810616458i;-0.029852749784765+0.382339810616458i,0.820818673768846+0i];
RR(:,:,2)=[0.820891385571210+0i,0.031198830973676-0.382170970063314i;0.031198830973676+0.382170970063314i,0.179108614428790+0i];



R_pure = zeros(d,d,k);
for i = 1:k
    [V,~] = eig(RR(:,:,i));
    psi = V(:,end);
    R_pure(:,:,i) = psi * psi';
end


ci_purif = compute_ci_purification(p, R_pure, K, n);
fprintf('   CI = %.14f\n\n', ci_purif);
