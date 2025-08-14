% comp_kraus.m
% This function computes the Kraus operators of the complementary channel for a given set of operators K
% and a specified output dimension.
% K ... cell array of Kraus operators
% dim ... vector specifying the output dimensions [d_out_B, d_out_E]



function [Kc] = comp_kraus(K,dim)
a = dim(1);
b = dim(2);

if (iscell(K))
    for k=1:length(K)
        kraus(:,:,k) = K{k};
    end
else
    kraus = K;
end     

for l=1:b
    Op(:,:) = kraus(l,:,:);
    Kc{l} = transpose(Op);
end
end