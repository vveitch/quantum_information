function [ mub_el ] = hw_mub( dim, m, r )
%returns the rth element of the mth HW covariant mutually unbiased basis 

if 0==mod(dim,2)
%for dim=power of 2, according to construction in arXiv:0704.1277v1 
%TODO: generalize construction to dimensions other than 4, allow other sets
%of HW covariant bases

    %hardcoded
    dim=4; 
    U=[
    1 1i 1i -1;
    1i 1 -1 1i;
    1 1i -1i 1;
    -1i -1 -1 1i;
    ]/2;

    %we'll assume indexing starts at 0
    r=mod(r,dim)+1;

    %rth comp basis projector
    mub_el=zeros(dim,1);
    mub_el(r)=1;

    mub_el=U^m*mub_el;
%    mub_el=mub_el'*mub_el;
elseif isprime(dim)
    %returns the rth element of the mth HW covariant mutually unbiased basis constructed according
    %to 'Spectra of phase point operators' paper. Dimension must be prime for
    %this to make any sense

    mub_el=zeros(dim);

    %the m=inf case from the paper (typical use: infinity=3)
    if m==dim
    for j=0:dim-1
        mub_el=mub_el+exp(-1i*2*pi/dim*r*j)*weyl_operators(dim,j,0);
    end

    % %HW operator cycles through mub
    % for k=1:dim-1
    %     mub(:,:,dim-k)=weyl_operators(dim,0,1)*mub(:,:,k+1)*weyl_operators(dim,0,1)';
    % end

    else
    for j=0:dim-1
        mub_el=mub_el+exp(-1i*2*pi/dim*r*j)*weyl_operators(dim,mod(m*j,dim),j);
    end

    % %HW operator cycles through mub
    % for k=1:dim-1
    %     mub(:,:,k+1)=weyl_operators(dim,1,0)*mub(:,:,k)*weyl_operators(dim,1,0)';
    % end
    end

    mub_el=mub_el/dim;
end