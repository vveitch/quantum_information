function [ weyl ] = prime_weyl_operators( dim, p, q )
%WEYL_OPERATORS Generates the prime dimension weyl operators as described in D. Gross'
%paper
%   Detailed explanation goes here
x=zeros(dim);
z=zeros(dim);
for l=1:dim-q        
    x(l+q,l)=1;
    z(l,l)=exp(p*(l-1)*2i*pi/dim);     
end
for l=dim-q+1:dim
    x(l-dim+q,l)=1;
    z(l,l)=exp(p*(l-1)*2i*pi/dim);     
end

weyl=exp(-p*q*(dim+1)*1i*pi/dim)*z*x;

end
