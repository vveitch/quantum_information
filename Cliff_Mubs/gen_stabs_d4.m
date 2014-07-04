%This code spits out 3 disjoint sets of MUBs composed of d(d+1)=20 vectors each
%these mubs are composed of stabilizer states (and use all 60 of them)

H=1/sqrt(2)*[1 1; 1 -1];
P=[1 0; 0 1i];
CNOT=[1 0 0 0;
    0 1 0 0;
    0 0 0 1;
    0 0 1 0];
SWAP=[1 0 0 0;
    0 0 1 0;
    0 1 0 0;
    0 0 0 1];
H_1=kron(H,eye(2));
P_2=kron(eye(2),P);

%these are the GHW mubs
MUBs=zeros(4,4,5);
for m=1:5
    for r=1:4
       MUBs(:,r,m)=hw_mub_vec(4,m-1,r-1);
    end
end

%a unitarily equivalent disjoint set
more_MUBs=zeros(4,4,5);
for m=1:5
    more_MUBs(:,:,m)=CNOT*SWAP*MUBs(:,:,m);
end

%another unitarily equivalent set, disjoint from the other 2 sets of mubs
even_more_MUBs=zeros(4,4,5);
for m=1:5
    even_more_MUBs(:,:,m)=(CNOT*SWAP)^2*MUBs(:,:,m);
end
