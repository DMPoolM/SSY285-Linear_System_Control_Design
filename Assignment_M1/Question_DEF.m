clc;
close all;
clear;

% give parameters value
R = 1;
K_E = 10^-1;
K_T = 10^-1;
J_1 = 10^-5;
J_2 = 4 * 10^-5;
B_f = 2 * 10^-3;
D_1 = 20;
D_2 = 2;

% define A matrics

A = [0,0,0,1,0;
    0,0,0,0,1;
    0,D_2/B_f,-D_2/B_f,0,0;
    -D_1/J_1,D_1/J_1,0,-(K_E*K_T)/(J_1*R),0;
    D_1/J_2,-(D_1+D_2)/J_2,D_2/J_2,0,0];
% calculate eigenvalues
eigenvalue = eig(A);

% Question D
% define B C D matrics 
B = [0,0;0,0;0,1/B_f;K_T/(J_1*R),0;0,0];
% for two cases of different choices of y(t)
C1 = [0,1,0,0,0;0,0,0,0,1];
C2 = [0,0,0,-K_E/R,0;0,D_2/B_f,-D_2/B_f,0,0];
D1 = [0,0;0,0];
D2 = [1/R,0;0,1/B_f];

% formulate G function
syms s;
I = eye(size(A));
%G1 = C1*(s*I-A)^-1 *B+D1;
G1 = ss(A, B, C1, D1);
%G2 = C2*(s*I-A)^-1 *B+D2;
G2 = ss(A, B, C2, D2);
G11 = tf(G1);
G22 = tf(G2);
% find zeros and poles

% case 1
for m = 1:1:2
    for n = 1:1:2
        [zero_a, pole_a, k_a] = tf2zp(G11.numerator{m,n}, G11.denominator{m,n})
    end
end
% case 2
for m = 1:1:2
    for n = 1:1:2
        [zero_b, pole_b, k_b] = tf2zp(G22.numerator{m,n}, G22.denominator{m,n})
    end
end
Ce = [0 1 0 0 0; 0 0 0 -0.1 0; 0 0 0 0 1; 0 1000 -1000 0 0];
De = [0 0; 1 0; 0 0; 0 500];

SS = ss(A,B,double(Ce),double(De))

s = tf('s')
I = eye(size(A));
G = tf((s*I-A)\B);

tv  = linspace(0,0.08,1000);
uv = [10*ones(size(tv));
      0*ones(1,numel(tv)/2) -0.1*ones(1,numel(tv)/2) ];

% step(G(2,1)*10)

% figure
Ce2 = [0 0 0 -0.1 0; 0 0 0 1 0; 0 0 0 0 1; 0 1000 -1000 0 0; 1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0];
De2 = [1 0; 0 0; 0 0; 0 500; 0 0; 0 0; 0 0];

SS = ss(A,B,double(Ce2),double(De2))

s = tf('s')
I = eye(size(A));
G = tf((s*I-A)\B);

tv  = linspace(0,0.08,1000);
uv = [10*ones(size(tv));
      0*ones(1,numel(tv)/2) -0.1*ones(1,numel(tv)/2) ];

% step(G(2,1)*10)

[ysol,tsol,xsol] = lsim(SS,uv,tv)

lw=2;
figure('Color','white')
subplot(9,1,1)
plot(tsol,uv(1,:), 'LineWidth', lw); grid on;
ylabel('v_a')
subplot(9,1,2)
plot(tsol,uv(2,:), 'LineWidth', lw); grid on;
ylabel('T_e')
subplot(9,1,3)
plot(tsol,ysol(:,1), 'LineWidth', lw); grid on;
ylabel('i_{a}')
subplot(9,1,4)
plot(tsol,ysol(:,2), 'LineWidth', lw); grid on;
ylabel('\omega{1}')
subplot(9,1,5)
plot(tsol,ysol(:,3), 'LineWidth', lw); grid on;
ylabel('\omega_{2}')
subplot(9,1,6)
plot(tsol,ysol(:,4), 'LineWidth', lw); grid on;
ylabel('\omega_{3}')
subplot(9,1,7)
plot(tsol,ysol(:,5), 'LineWidth', lw); grid on;
ylabel('\phi_{1}')
subplot(9,1,8)
plot(tsol,ysol(:,6), 'LineWidth', lw); grid on;
ylabel('\phi_{2}')
subplot(9,1,9)
plot(tsol,ysol(:,7), 'LineWidth', lw); grid on;
ylabel('\phi_{3}')


xlabel('time [s]')
set(gca,'LooseInset',get(gca,'TightInset'))
saveas(gcf, fullfile(pwd,'images/1e_outputs'),'epsc')