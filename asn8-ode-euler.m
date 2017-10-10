clear all 
close all
clc


% initiate the problem 

% ivp=[phi1; phi1_dot; phi2; phi2_dot; g; m1; m2; l1; l2] - vector of
% parameters and solution of ODE

% tend - termination time of simulations

x0 = [pi/6; 0; pi/6; 0; 9.8; 2; 1; 2; 1];
tend = 10;
average_runningTimes = zeros(1, 3);
runningTime = zeros(1, 10);
numTrials = 20;


%% ODE45 solution

for k = 1 : numTrials
    
    tic;    % start timer
    [t45, sol_ode45] = ode45( @double_pendulum_ODE, [0 tend], x0 );
    runningTime( 1 ) = toc;   % save running time

end

average_runningTimes( 1 ) = mean( runningTime );


%% Forward Euler solution

h = [10.^(-1) 10.^(-3) ];%10.^(-4)];

for j = 1 : numel(h)
    
    for k = 1 : numTrials
        
        sol_Euler(:,1) = x0;
        time = 0:h(j):tend;

        % Put code here to solve using forward Euler method.
        tic;
        for i = 1:( numel(time) - 1 )
            sol_Euler(:, i + 1) = sol_Euler(:, i) + h(j) * double_pendulum_ODE( transpose(time), sol_Euler(:, i) );
        end

        runningTime( k ) = toc;
    end

    if j == 1
        solution_Euler_h_1 =  sol_Euler;
    elseif j == 2
        solution_Euler_h_2 =  sol_Euler;
%     else % if j == 3
%         solution_Euler_h_3 =  sol_Euler;
    end
    
    average_runningTimes( j + 1 ) = mean( runningTime );
    
end


%% plot solutions

figNum = 1;  % start with figure 1

for k = 1:2
    figure( figNum );
    
    if k == 1
        plot(t45, sol_ode45(:,1), 'r', ...
            0:h(1):tend, solution_Euler_h_1(1, 1:end ), 'g', ...
            0:h(2):tend, solution_Euler_h_2(1, 1:end ), 'b' );%, ...
%             0:h(3):tend, solution_Euler_h_3(1, 1:end ), 'c');    
    else % if k == 2
        plot(t45, sol_ode45(:,3), 'r' , ...
            0:h(1):tend, solution_Euler_h_1(3, 1:end ), 'g', ...
            0:h(2):tend, solution_Euler_h_2(3, 1:end ), 'b' );%, ...
%             0:h(3):tend, solution_Euler_h_3(3, 1:end ), 'c');
    end
    
    axis([0 tend -1 1]) % limit y-range to -1 < degree < 1
    title( {['Figure ' num2str( figNum ) ':  Approximated \phi_' num2str( k ) 's to Double Pendulum Problem'], ...
            ['when  0 <= t <= ' num2str( tend )]}, 'fontsize', 10 )
    xlabel('Time (t)')
    ylabel('angle')
    legend('using ODE45', ['Euler h = ' num2str( h(1) )], ['Euler h = ' num2str( h(2) )], 'Location', 'northwest', 'Orientation','horizontal');%, ['Euler h = ' num2str( h(3) )])
    figNum = figNum + 1;
    
end


%% comparing efficiency
figure( figNum );
bar( average_runningTimes(1:3) );
set(gca,'xticklabel',{'ODE45',['Euler h = ' num2str( h(1) )], ['Euler h = ' num2str( h(2) )]});%, ['Euler h = ' num2str( h(3) )]});
title( {['Figure ' num2str( figNum ) ':  Average computational time required to solve '], ...
        ['the Double Pendulum Problem (after ', num2str( numTrials ), ' iterations)' ]}, 'fontsize', 10 )
xlabel('Method')
ylabel('Running Time (seconds)')



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
