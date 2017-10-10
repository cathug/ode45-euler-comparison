clear all 
close all
clc


% initiate the problem 

% ivp=[phi1; phi1_dot; phi2; phi2_dot; g; m1; m2; l1; l2] - vector of
% parameters and solution of ODE

% tend - termination time of simulations

x0 = [pi/6; 0; pi/6; 0; 9.8; 2; 1; 2; 1];
tend = 200;
fignum = 4;


%% ODE45 solution

[t45, sol_ode45] = ode45( @double_pendulum_ODE, [0 tend], x0 );


%% Forward Euler solution

h = 10.^(-3);

sol_Euler(:,1) = x0;
time = 0:h:tend;

% Put code here to solve using forward Euler method.
for i = 1:( numel(time) - 1 )
    sol_Euler(:, i + 1) = sol_Euler(:, i) + h * double_pendulum_ODE( transpose(time), sol_Euler(:, i) );
end



% compute total energy of solutions 

m1 = x0(6); m2 = x0(7); l1= x0(8); l2=x0(9); g = x0(5);
phi_1_ode45 = sol_ode45(:,1); dphi_1_dt_ode45 = sol_ode45(:,2);
phi_2_ode45 = sol_ode45(:,3); dphi_2_dt_ode45 = sol_ode45(:,4);
phi_1_euler = sol_Euler(1, 1:end ); dphi_1_dt_euler = sol_Euler(2, 1:end );
phi_2_euler = sol_Euler(3, 1:end ); dphi_2_dt_euler = sol_Euler(1, 1:end );


% Here you need to correctly compute the potential and kinetic energy with
% your solution from ode45

P_ode45 = -g * ( ( m1 + m2 ) * l1 * cos( phi_1_ode45 ) + m2 * l2 * cos ( phi_2_ode45 ) );
prod_1 = dphi_1_dt_ode45 * l1;
prod_2 = dphi_2_dt_ode45 * l2;
K_ode45 = 0.5 * ( (m1 + m2) * prod_1.^2 + m2 * prod_2.^2 ) + m2 * prod_1 .* prod_2 .* ...
    ( cos( phi_1_ode45 .* phi_2_ode45 ) + sin( phi_1_ode45 .* phi_2_ode45 ) ) ;

figure( fignum )
plot(t45, P_ode45 + K_ode45)
xlabel('time')
ylabel('Total energy')
title( {['Figure ' num2str( fignum ) ':  Estimated total energy of the system using ode45'], ...
    ['when 0 <= t <= ' num2str( tend )]}, 'fontsize', 10 )
fignum = fignum + 1;

% Here you need to correctly compute the potential and kinetic energy with
% your solution from forward Euler.

P_Euler =  -g * ( ( m1 + m2 ) * l1 * cos( phi_1_euler ) + m2 * l2 * cos ( phi_2_euler ) );
prod_1 = dphi_1_dt_euler * l1;
prod_2 = dphi_2_dt_euler * l2;
K_Euler = 0.5 * ( (m1 + m2) * prod_1.^2 + m2 * prod_2.^2 ) + m2 * prod_1 .* prod_2 .* ...
    ( cos( phi_1_euler .* phi_2_euler ) +  sin( phi_1_euler .* phi_2_euler ) );

figure( fignum )
plot(time, P_Euler + K_Euler)
xlabel('time')
ylabel('Total energy')
title( {['Figure ' num2str( fignum ) ':  Estimated total energy of the system using Eulers method'], ...
    ['when 0 <= t <= ' num2str( tend )]}, 'fontsize', 10 )





% this function defines the double pendulum ODE 
function xdot = double_pendulum_ODE(t,x)

%
%   author:  Alexander Erlich (alexander.erlich@gmail.com)
%
%   parameters:
%
%   t       Column vector of time points 
%   xdot    Solution array. Each row in xdot corresponds to the solution at a
%           time returned in the corresponding row of t.
%

    g=x(5); m1=x(6); m2=x(7); l1=x(8); l2=x(9);

    xdot=zeros(9,1);

    xdot(1)=x(2);

    xdot(2)=-((g*(2*m1+m2)*sin(x(1))+m2*(g*sin(x(1)-2*x(3))+2*(l2*x(4)^2+...
        l1*x(2)^2*cos(x(1)-x(3)))*sin(x(1)-x(3))))/...
        (2*l1*(m1+m2-m2*cos(x(1)-x(3))^2)));

    xdot(3)=x(4);

    xdot(4)=(((m1+m2)*(l1*x(2)^2+g*cos(x(1)))+l2*m2*x(4)^2*cos(x(1)-x(3)))*...
        sin(x(1)-x(3)))/(l2*(m1+m2-m2*cos(x(1)-x(3))^2));

end
