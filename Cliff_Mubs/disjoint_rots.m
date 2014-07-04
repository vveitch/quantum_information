function [ disj_rots ] = disjoint_rots( n )
%DISJOINT_ROTS this returns a set of unitaries {U} such that the mubs
%{{R^i}*U_j*comp_basis} are all disjoint
%n must be at least 2

%basically these U_j are {R_n}{R_n-1}..{R_1}

prim_rot=prim_rots(n-1);

if n~=2
    %recursive bit
    dr_small=disjoint_rots(n-1);
    small=size(dr_small,3);
    
    disj_rots=zeros(2^n,2^n,(2^(n-1)+1)*small);
    
    for i=1:2^(n-1)+1
        for j=1:small
            %note the identity on the first qubit
            disj_rots(:,:,small*(i-1)+j)=kron(eye(2),prim_rot^(i-1)*dr_small(:,:,j));
        end
    end
else
    %base case
    disj_rots=zeros(4,4,3);
    for i=1:3
        disj_rots(:,:,i)=kron(eye(2),prim_rot^i);
    end
end

end
