function [ meas ] = tensor_operators( cell_Ops )
%TENSOR_OPERATORS tensor every gorram thing
%   Take the tensor product of everything with everything else in order to
%   create a tomographically complete basis for a composite Hilbert space
%   cell_Ops should contain an array of povms, each of which is an array of
%   projectors; should look like (subsystem_index, operator_index, matrix
%   indices [l, m])


sub_ops_a=cell2mat(cell_Ops(1));
dim_a=length(sub_ops_a(1,1,:));
if length(cell_Ops)>2
    %figure out the final size of what I'm spitting out for preallocation
    meas_size=dim_a;
    for j=2:length(cell_Ops)
        meas_size=meas_size*length(cell_Ops{j});
    end
    meas=zeros(fix(sqrt(meas_size)),fix(sqrt(meas_size)),meas_size);
    %the tensor thing is easy enough to figure out as long as you've only
    %ever got 2 of this freaking things going at once
    sub_ops_b=tensor_operators(cell_Ops(2:end));
    dim_b=length(sub_ops_b(1,1,:));
    for j=1:dim_a
        for k=1:dim_b
            meas(:,:,(j-1)*dim_b+k)=kron(sub_ops_a(:,:,j),sub_ops_b(:,:,k));
        end
    end
    %if you've only got 2 povms to tensor then this is your guy
elseif length(cell_Ops)==2
    sub_ops_b=cell2mat(cell_Ops(2));
    dim_b=length(sub_ops_b(1,1,:));
    meas_size=dim_a*dim_b;
    meas=zeros(fix(sqrt(meas_size)),fix(sqrt(meas_size)),meas_size);
    for j=1:dim_a
        for k=1:dim_b
            meas(:,:,(j-1)*dim_b+k)=kron(sub_ops_a(:,:,j),sub_ops_b(:,:,k));
        end
    end
end
end
