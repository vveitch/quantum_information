function rot = prim_rots(qubits)
%this returns a unitary matrix which cycles through full sets of MUBs

%TODO: the prim rot doesn't have the same generator decomposition in all
%dimensions because god hates me. Gotta find a way around that

%we exploit the mathematics of galois fields ala
%http://arxiv.org/PS_cache/arxiv/pdf/0704/0704.1277v1.pdf to do this
%R=L_1*L_2*L_3 - we want the unitary rep of this
%representations of L_1, L_2 and L_3 are given in Gibbons, Hoffman, Wooters

if qubits==1
    rot=[0.5000 + 0.5000i   0.5000 + 0.5000i;
        0.5000 - 0.5000i  -0.5000 + 0.5000i];
elseif qubits==2
    rot=[0.5000                  0 + 0.5000i        0 - 0.5000i   0.5000    ;
        0 + 0.5000i   0.5000            -0.5000                  0 + 0.5000i;
        0.5000                  0 - 0.5000i        0 - 0.5000i  -0.5000          ;
        0 - 0.5000i   0.5000             0.5000                  0 + 0.5000i];
elseif qubits==3
    rot=[
        0.3536   -0.3536   -0.3536    0.3536   -0.3536    0.3536    0.3536   -0.3536;
        0 - 0.3536i        0 - 0.3536i        0 - 0.3536i        0 - 0.3536i         0 - 0.3536i        0 - 0.3536i        0 - 0.3536i        0 - 0.3536i;
        -0.3536   -0.3536    0.3536    0.3536   -0.3536   -0.3536    0.3536    0.3536;
        0 + 0.3536i        0 - 0.3536i        0 + 0.3536i        0 - 0.3536i         0 - 0.3536i        0 + 0.3536i        0 - 0.3536i        0 + 0.3536i;
        -0.3536    0.3536   -0.3536    0.3536   -0.3536    0.3536   -0.3536    0.3536;
        0 + 0.3536i        0 + 0.3536i        0 - 0.3536i        0 - 0.3536i         0 - 0.3536i        0 - 0.3536i        0 + 0.3536i        0 + 0.3536i;
        -0.3536   -0.3536   -0.3536   -0.3536    0.3536    0.3536    0.3536    0.3536;
        0 + 0.3536i        0 - 0.3536i        0 - 0.3536i        0 + 0.3536i         0 + 0.3536i        0 - 0.3536i        0 - 0.3536i        0 + 0.3536i];
   
else
    CNOT=[1 0 0 0;
        0 1 0 0;
        0 0 0 1;
        0 0 1 0];
    SWAP=[1 0 0 0;
        0 0 1 0;
        0 1 0 0;
        0 0 0 1];
    l_1=[1 0; 0 1i];
    l_2=1/sqrt(2)*[1 1i; 1i 1];
    
    rot=kron(SWAP,eye(2^(qubits-2)));
    for i=1:qubits-2
        rot=kron(eye(2^i),kron(SWAP,eye(2^(qubits-2-i))))*rot;
    end
    
    %I'm kind of fucking the pooch on the decomposition in the GHW paper, no compelling
    %reason to believe this is generally an order d-1 operator
    %Actually, since *nothing* works for 3 qubits I'm inclined to suspect
    %the GHW paper is wrong
    rot=kron(CNOT,eye(2^(qubits-2)))*rot;
    %rot=blkdiag(eye(4),kron(eye(2),[0 1; 1 0]))*rot;
    
    %this creates the reps of the L_1 and L_2 matrices
    L_1=l_1;
    L_2=l_2;
    for i=2:qubits
        L_1=kron(l_1,L_1);
        L_2=kron(l_2,L_2);
    end
    %this decomposition is total numerical black magic
    %but I think L_2*L_1*L_3*L_2*L_1 would work if I had the earlier decomp
    %right
    rot=rot^(2^(qubits-1)+1)*L_2*L_1;
end
end