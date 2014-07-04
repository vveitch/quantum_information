function [ rho ] = recon( rep )
%RHO_NEW This takes the cov rep of a state and outputs the associated
%density matrix
%dimension is assumed prime for now

dim=size(rep,1);
rho=zeros(dim);
ppos=cov_phase_point_ops(dim,1);

for i=1:dim
    for j=1:dim
        rho=rho+rep(i,j)*ppos(:,:,i,j);
    end
end

end

