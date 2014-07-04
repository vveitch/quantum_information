function [ phase_point ] = phase_point_operators_take_two( dim, q, p )
%returns the phase point operators for wigner functions according to the
%formula given in D. Gross' paper.
phase_point=zeros(dim);

for j=0:dim-1
    for k=0:dim-1
   
       phase_point=phase_point+exp(1i*2*pi/dim*(q*k-p*j)) ...
           *weyl_operators(dim,j,k)';
    end
end
phase_point=phase_point/dim;

