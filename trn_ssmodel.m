function [A,B,C,D ] = trn_ssmodel(p,nstates,ninputs,alpha,cz)
% Construct A, B, C and D matrices for the zone

A = zeros(nstates);
B = zeros(nstates,ninputs);
C = zeros(1,nstates);
D = zeros(1,ninputs);

ae = alpha(1);
ag = alpha(2);


A(1,1) = (-p(1)-p(2))/p(8);
A(1,2) = p(2)/p(8);

A(2,1) = p(2)/p(9);
A(2,2) = (-p(2)-p(3))/p(9);
A(2,5) = p(3)/p(9);

A(3,3) = (-p(4)-p(5))/p(10);
A(3,4) = p(5)/p(10);

A(4,3) = p(5)/p(11);
A(4,4) = (-p(5)-p(6))/p(11);
A(4,5) = p(6)/p(11);

A(5,2) = p(3)/cz;
A(5,4) = p(6)/cz;
A(5,5) = (-p(3)-p(6)-p(7))/cz;


B(1,1) = p(1)/p(8);
B(1,3) = p(3)/p(8);
B(2,4) = ae/p(9);
B(2,5) = ae/p(9);
B(3,2) = p(2)/p(10);
B(4,4) = ag/p(11);
B(4,5) = ag/p(11);
B(5,1) = p(7)/cz;
B(5,6) = 1/cz;
B(5,7)= -1/cz;

C(1,5) = 1;

end

