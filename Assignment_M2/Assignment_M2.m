clc;
close all;
clear;

% give parameters value
%R = 1;
K_E = 10^-1;
K_T = 10^-1;
J_1 = 10^-5;
J_2 = 4 * 10^-5;
B_f = 2 * 10^-3;
%D_1 = 20;
D_2 = 2;
syms R D_1
% define A matrics

A = [0,0,0,1,0;
    0,0,0,0,1;
    0,D_2/B_f,-D_2/B_f,0,0;
    -D_1/J_1,D_1/J_1,0,-(K_E*K_T)/(J_1*R),0;
    D_1/J_2,-(D_1+D_2)/J_2,D_2/J_2,0,0];
% calculate eigenvalues
eigenvalue = eig(A);

% Question C
% define B C D matrics 
B = [0,0;0,0;0,1/B_f;K_T/(J_1*R),0;0,0];
% for two cases of different choices of y(t)
C1 = [0,1,0,0,0;0,0,0,0,1];
C2 = [0,0,0,-K_E/R,0;0,D_2/B_f,-D_2/B_f,0,0];
D1 = [0,0;0,0];
D2 = [1/R,0;0,1/B_f];
W_c = [B,A*B,A^2*B,A^3*B,A^4*B]
r_c = rref(W_c)
W_o1 = [C1; C1*A; C1*A^2; C1*A^3; C1*A^4]
W_o2 = [C2; C2*A; C2*A^2; C2*A^3; C2*A^4]
r_o1 = rref(W_o1)
r_o2 = rref(W_o2)
simplify (r_c)
simplify (r_o1)
simplify (r_o2)
latex(r_c)

I = eye(size(A));
PBH1 = [eigenvalue(1)*I - A; C1]
PBH2 = [eigenvalue(1)*I - A; C2]
for i = 1:5
    PBH3 = [eigenvalue(i)*I - A, B];
    a = rank(PBH3)
end
rref(PBH1)
rref(PBH2)
% latex(rref(PBH2))
% latex (eigenvalue)