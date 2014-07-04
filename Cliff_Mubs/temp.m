% %some useful operators
% H=1/sqrt(2)*[1 1; 1 -1];
% P=[1 0; 0 1i];
% CNOT=[1 0 0 0;
%     0 1 0 0;
%     0 0 0 1;
%     0 0 1 0];
% SWAP=[1 0 0 0;
%     0 0 1 0;
%     0 1 0 0;
%     0 0 0 1];
% H_1=kron(H,eye(2));
% P_2=kron(eye(2),P);
% %paulis
% Z=[1 0; 0 -1];
% X=[0 1; 1 0];
%
% p_1=kron(Z,Z);
% p_2=kron(eye(2),Z);
% p_3=kron(Z,eye(2));
% A = [p_1 p_2 p_3];
% [V, D] = joint_diag(A,1.0e-8);


stabs=d4;
stabs=round(100000*stabs)/100000;

%indices of the special MUB representation
ind = [1 2 3 5 6 7 9 10 11 13 14 15 17 18 19];
specialMUB = stabs(:,ind);
%represent them all in the special MUB
for state = 1:length(stabs);
    for component = 1:15
        MUBrep(state,component) = abs(stabs(:,state)'*specialMUB(:,component))^2;
    end
end
foo = round(100*MUBrep)/100;
unique(foo,'rows')
