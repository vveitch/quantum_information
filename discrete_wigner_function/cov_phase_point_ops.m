function [ ppos ] = cov_phase_point_ops( base_dim, exponent )
%COV_PHASE_POINT_OPS this spits out the phase point operators which have
%non-negative rep for all stabilizer states
%These are covariant under the action of the Clifford group for d odd
%For d a power of 2 they are covariant up to complex conjugation of the
%phase point operators tensor product wise


ppos=zeros(base_dim^exponent,base_dim^exponent,base_dim^exponent,base_dim^exponent);

weyl_vec=zeros(2, exponent);

for i=0:base_dim^(2*exponent)-1
    %this just produces every vector of the form [a b][c d]..[l k] where
    %the 2n entries in each vector are taken from {0..base_dim}
    blech=i;
    for j=1:exponent
        weyl_vec(1,j)=mod(blech,base_dim);
        blech=(blech-weyl_vec(1,j))/base_dim;
        weyl_vec(2,j)=mod(blech, base_dim);
        blech=(blech-weyl_vec(2,j))/base_dim;
    end
    %the 0,0 ppo is the sum of the HW operators
    ppos(:,:,1,1)=ppos(:,:,1,1)+weyl_operators(base_dim,exponent,weyl_vec);
end

ppos(:,:,1,1)=ppos(:,:,1,1)/(base_dim^exponent);

%exploit the HW covariance
for i=0:base_dim^(2*exponent)-1
    %this just produces every vector of the form [a b][c d]..[l k] where
    %the 2n entries in each vector are taken from {0..base_dim}
    blech=i;
    p=0;
    q=0;
    for j=1:exponent
        weyl_vec(1,j)=mod(blech,base_dim);
        p=p+weyl_vec(1,j)*base_dim^(j-1);
        blech=(blech-weyl_vec(1,j))/base_dim;
        weyl_vec(2,j)=mod(blech, base_dim);
        q=q+weyl_vec(2,j)*base_dim^(j-1);
        blech=(blech-weyl_vec(2,j))/base_dim;
    end
    %the 0,0 ppo is the sum of the HW operators
    ppos(:,:,p+1,q+1)=weyl_operators(base_dim,exponent,weyl_vec)'*ppos(:,:,1,1)*weyl_operators(base_dim,exponent,weyl_vec);
end
