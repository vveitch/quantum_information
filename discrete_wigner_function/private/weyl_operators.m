function [ weyl ] = weyl_operators( base_dim, exponent, op_vec)
%weyl_operators this returns the power of prime weyl operator corresponding to the op_vec
%(the 2n vector in f_base dim in the symplectic rep)
%   base_dim is the prime number which is the dimension of the base system
%   (ie. qubit, qutrit, qupit)
%   exponent is the number of qupits
%   op_vec is the 2n vector in f_base dim in the symplectic rep, it should
%   be fed in as a matrix of the dimensionality [2 exponent] so that
%   each column corresponds to a term in the direct sum

weyl=prime_weyl_operators(base_dim, op_vec(1,1), op_vec(2,1));
for i=2:exponent
    weyl=kron(weyl,prime_weyl_operators(base_dim,op_vec(1,i),op_vec(2,i)));
end


end

