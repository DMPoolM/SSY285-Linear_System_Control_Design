clear
%% 3-DOF helicopter nonlinear model parameters %%

Mf = 0.71; %[kg] Mass of front propeller assembly
Mb = 0.71; %[kg] Mass of back propeller assembly
Mc = 1.69; %[kg] Mass of the counterwight
Ld = 0.05; %[m] Length of pendulum for elevention axis
Lc = 0.44; %[m] Distance from pivot point to counterweight
La = 0.62; %[m] Distance from pivot point to helicopter body
Le = 0.02; %[m] Length of pendulum for pitch axis
Lh = 0.18; %[m] Distance from pitch axis to motor
g = 9.81; %[m/s^2] Gravitational acceleration
Km = 0.12; %[N/V] Propeller force-thrust constant
Jeps = 0.86; %[kg*m^2] Moment of inertia level axis
Jro = 0.044; %[kg*m^2] Moment of inertia pitch axis
Jlam = 0.82; %[kg*m^2] Moment of inertia travel axis
Neps = 0.001; %[kg*m^2/s] Viscous friction level axis
Nro = 0.001; %[kg*m^2/s] Viscous friction pitch axis
Nlam = 0.005; %[kg*m^2/s] Viscous friction travel axis

% Coefficients p
p1 = (-(Mf+Mb)*g*La+Mc*g*Lc)/Jeps;
p2 = (-(Mf+Mb)*g*(Ld+Le)+Mc*g*Ld)/Jeps;
p3 = -Neps/Jeps;
p4 = Km*La/Jeps;
p5 = (-Mf+Mb)*g*Lh/Jro;
p6 = -(Mf+Mb)*g*Le/Jro;
p7 = -Nro/Jro;
p8 = Km*Lh/Jro;
p9 = -Nlam/Jlam;
p10 = -Km*La/Jlam;

% initial conditions
eps0 = 0;
ro0 = 0;
lam0 = 0;
eps0dot = 0;
ro0dot = 0;
lam0dot = 0;
x0 = [eps0; ro0; lam0; eps0dot; ro0dot; lam0dot]
%Linearize
syms Vf Vb
u = [Vf; Vb];
fx0 = [eps0dot; ro0dot; lam0dot; p1*cos(eps0)+ p2 * sin(eps0)+ p3* eps0dot; p5*cos(ro0)+ p6 * sin(ro0)+ p7* ro0dot;p9* lam0dot];
gx0 = [0 0; 0 0; 0 0; p4* cos(lam0) p4* cos(lam0); p8 -p8; p10* sin(ro0) p10* sin(ro0)];
[Vf0, Vb0]=solve(fx0 + gx0*u == 0);
Vf0 = double(Vf0);
Vb0 = double(Vb0);
u0 = [Vf0;Vb0];
A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     -0.169963953488372 0 0 -0.00116279069767442 0 0;
     0 -6.33190909090909 0 0 -0.0227272727272727 0;
     0 -1.63659512195122 0 0 0 -0.00609756097560976]
B = [0 0;
     0 0;
     0 0;
     0.0865116279069767 0.0865116279069767;
     0.490909090909091 -0.490909090909091;
     0 0]
C = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0]
D = [0 0; 0 0; 0 0]

% Vf0 = 1;
% Vb0 = 0;

sim_out = sim('linearized_model', 'StopTime', '10');
plot_results(sim_out.sim_x_lin, sim_out.sim_x_nonlin, sim_out.sim_u)
sys = ss(A, B, C, D);


%p = [-1 -2 -3 -4 -5 -6];
%K = place(A, B, p)
% LQ
% Qx = diag([1 1 1 1 1 1]);
% Qu = eye(length(u));
% [K, S, e] = lqr(A, B, Qx, Qu)
% Kr = pinv(C*inv(B*K-A)*B);
% t_R = [1; 3];
% R = [50; 50]/180*pi;
% sim_out = sim('LQ', 'StopTime', '10');
% plot_results(sim_out.sim_x_lin, sim_out.sim_x_nonlin, sim_out.sim_u)

% LQI
Ci = [1 0 0 0 0 0;
     0 0 1 0 0 0];
A_I = [A Ci'*0 ;
      -Ci Ci*Ci'*0];
B_I = [B; Ci*Ci'*0];

Qx = diag([1 1 1 1 1 1 1 1]);
Qu = eye(length(u))/100;

K_lqr = lqr(A_I, B_I, Qx, Qu);

K   = K_lqr(:,1:6);
K_I = K_lqr(:,6+1:end);
Kr = pinv(C*inv(B*K-A)*B);

t_R = [0; 2];
R = [0; 360]/180*pi;
sim_out = sim('LQI', 'StopTime', '20');
plot_results(sim_out.sim_x_lin, sim_out.sim_x_nonlin, sim_out.sim_u)

function plot_results(sim_x_lin, sim_x_nonlin, sim_u)
    plot_sim = @(sim_out,index,style) plot(sim_out.Time, sim_out.Data(:,index)/pi*180, style, 'LineWidth',2);
    ylabels_X = {'level eps','pitch ro','travel lam','eps_{dot}','ro_{dot}','lam_{dot}'};
    plotorder = [1,3,5,2,4,6];
    figure('Color','white')
    for i=1:size(sim_x_lin.Data,2)
        subplot(4,2,plotorder(i)); hold on;
        plot_sim(sim_x_lin,i,'-');
        plot_sim(sim_x_nonlin,i,'--');
        xlabel('Time [s]')
        ylabel(ylabels_X{i})
    end
    legend({'Linear Model','Nonlinear Model'})
    
    ylabels_u = ({'Vf','Vb'});
    subplot(4,2,7)
    plot(sim_u.Time, sim_u.Data(:,1));  
    xlabel('Time [s]')
    ylabel(ylabels_u{1})
    subplot(4,2,8)
    plot(sim_u.Time, sim_u.Data(:,2));
    xlabel('Time [s]')
    ylabel(ylabels_u{2})
end