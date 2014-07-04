function MUB = d4
%this code spits out 9 full sets of MUBs (20 vectors each) composed of
%stabilizer states.
%these are divided into disjoint sets of 3
%See my notes "On the Structure of the Clifford MUBs"

%some useful operators
H=1/sqrt(2)*[1 1; 1 -1];
P=[1 0; 0 1i];
CNOT=[1 0 0 0;
    0 1 0 0;
    0 0 0 1;
    0 0 1 0];
H_1=kron(H,eye(2));
P_2=kron(eye(2),P);
R_1=1i*kron(eye(2),H*P*H*P);
%symplectic primitive rotation cycles through basis sets in mubs
U=[
    1 1i 1i -1;
    1i 1 -1 1i;
    1 1i -1i 1;
    -1i -1 -1 1i;
    ]/2;
%one of the generators for SL(2,F_4) from GHW paper; used to get alt rot
L_1=kron([1 0; 0 1i],[1 0; 0 1i]);

%%%GHW Mubs%%%
V=eye(4);
%these are the GHW mubs
MUBs=zeros(4,4,15);
for m=1:5
    for r=1:4

%the ops here are aribtrary members of the Clifford group which are not a multiple
%of the rotation operator up to permutation of the computational basis
        MUBs(:,r,m)=U^(m-1)*V(:,r);
        MUBs(:,r,m+5)=U^(m-1)*R_1*V(:,r);
        MUBs(:,r,m+10)=U^(m-1)*R_1^2*V(:,r);
    end
end
MUB = reshape(MUBs,4,60);

%the other choice of primitive rotation
U=L_1*U;

for m=1:5
    for r=1:4
        MUBs(:,r,m)=U^(m-1)*V(:,r);
        MUBs(:,r,m+5)=U^(m-1)*H_1*V(:,r);
        %notice that to get a 'good' clifford here CNOT*H_1 doesn't work,
        %because the primitive rotation is different
        MUBs(:,r,m+10)=U^(m-1)*P_2*CNOT*H_1*P_2*CNOT*H_1*V(:,r);
    end
end
MUB=[MUB reshape(MUBs,4,60)];

%standardize global phase
for m=1:length(MUB)
    for j=1:4
        if MUB(j,m)~=0
            MUB(:,m)=abs(MUB(j,m))/MUB(j,m)*MUB(:,m);
            break;
        end
    end
end

%%% Following, GHW we have freedom to choose the factor of the dual basis
%%% which changes the rep of the max commuting subgroups. These correspond
%%% to the 2 other choices. It turns out (numerically, and I think I can
%%% prove it) that this operation just cycles the bases in the MUB. Ie it
%%% doesn't affect the associated polytope at all.
% %%2nd way to partition stabilizers into MUBs
% %maximal commuting set of paulis
% p_1=kron(Z,X*Z);
% p_2=kron(X,Z);
% p_3=kron(X*Z,X);
%
% A = [p_1 p_2 p_3];
% [V, D] = joint_diag(A,1.0e-8);

% %%3rd way to partition stabilizers into MUBs
% %maximal commuting set of paulis
% p_1=kron(eye(2),X*Z);
% p_2=kron(X*Z,eye(2));
% p_3=kron(X*Z,X*Z);
%
% A = [p_1 p_2 p_3];
% [V, D] = joint_diag(A,1.0e-8);
