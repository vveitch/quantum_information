function [ weyl ] = weyl_operators( dim, q, p )
%WEYL_OPERATORS Generates the weyl operators as described in D. Gross'
%paper
%   Detailed explanation goes here
basis1=zeros(dim,1);basis2=zeros(dim,1);
weyl=zeros(dim);
for l=1:(dim)
    
    basis1(l)=1;
    basis2(mod(l-1+q,dim)+1)=1;
    
    weyl=weyl+exp(p*(l-1)*2*pi*1i/dim)*basis2*basis1';
    
    basis1(l)=0;
    basis2(mod(l-1+q,dim)+1)=0;
    
end

weyl=exp(-(dim-1)*p*q*1i*pi/dim)*weyl;

end
