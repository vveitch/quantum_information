function [ ppf ] = phase_point_frame( dim )

%PHASE_POINT_BASIS the phase point operators for a given dimension

%this code is really ineffecient, but it doesn't have to be run much

prime_decomp=factor(dim);

for h=1:length(prime_decomp)%composite dimensions made up of tensor products
    sub_dim=prime_decomp(h);
    sub_ppf=zeros(sub_dim,sub_dim,sub_dim, sub_dim);
    for j=0:sub_dim-1
        for k=0:sub_dim-1
            sub_ppf(:,:,j+1,k+1)=phase_point_operators_take_two(sub_dim,j,k)/sub_dim; %Chris changed this Dec05/2010
        end
    end
    C{h}=sub_ppf;
end

if length(C)==1 
    ppf=cell2mat(C(1));
else
    ppf=tensor_operators(C);
end


end

