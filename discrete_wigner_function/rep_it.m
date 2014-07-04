function [ wig_rho ] = rep_it( rho)
%REP_IT this will spit out the discrete Wigner representation of rho
% only works for power of prime dimension

d=size(rho,1);

wig_rho=zeros(d,d);

factors=factor(d);

ppos=cov_phase_point_ops(factors(1),size(factors,2))/d;

for i=1:d
    for j=1:d
        wig_rho(i,j)=real(trace(rho*ppos(:,:,i,j))); %this value is real anyways, but real() cleans up the output
%         wig_rho(i,j)=trace(rho*phase_point_operators_take_two(d,i-1,j-1));
    end
end


end

