% Rimless wheel simulation
clear all; close all;

% Parameters
params.g = 9.81;      % gravitational acceleration (m/s^2)
params.l = 1;         % leg length (m)
params.gamma = 0.08;   % slope angle (radians)
params.n = 8;         % number of spokes
params.alpha = 2*pi/params.n; % angle between spokes


% theta0 = 0; 
% omega0 = 4*sqrt((2*params.g/params.l)*(1-cos(params.gamma-params.alpha)));              % initial angular velocity
theta0 = -0.3;
omega0 = 1.5;

x0 = 0;                  % initial contact point x position
y0 = [theta0; omega0; x0];  % state vector
tspan = [0 4];          % simulation timespan

% Collision even stops the ode and gives back the time and trajectory of
% the segment
options = odeset('Events', @(t,y)collisionEvent(t,y,params), ...
                'RelTol', 1e-12, 'AbsTol', 1e-12, ...
                'MaxStep', 0.1);

% Initialize storage for complete trajectory
T = [];
Y = [];
t = tspan(1);
y = y0;

% Main simulation loop
while t(end) < tspan(2)
    % Simulate until next collision
    [t_temp, y_temp, te, ye, ie] = ode45(@(t,y)dynamics(t,y,params), ...
        [t(end) tspan(2)], y, options);
    
    % Store results
    if ~isempty(t_temp)
        if isempty(T)
            T = t_temp;
            Y = y_temp;
        else
            T = [T; t_temp(2:end)];
            Y = [Y; y_temp(2:end,:)];
        end
    end
    
    if isempty(ie)
        break
    end
    
    % Apply collision map (The starting point of the next rotation)
    y_plus = collisionMap(ye', params);
    t = te;
    y = y_plus;
end

% Animation
figure('Position', [100 100 800 600]);
for i = 1:10:length(T)
    clf;
    hold on;
    grid on;
    

    theta = Y(i,1);
    x_contact = Y(i,3);
    
    % Calculate ground height at contact point
    y_contact = x_contact * tan(params.gamma);
    
    % Calculate hub position relative to contact point
    x_hub = x_contact - params.l*sin(theta);
    y_hub = y_contact + params.l*cos(theta);
    
    % Draw slope
    x_ground = [x_hub-3 x_hub+3];
    y_ground = x_ground * tan(params.gamma);
    plot(x_ground, y_ground, 'k', 'LineWidth', 2);
    
    % Draw spokes
    for j = 1:params.n
        spoke_angle = theta + (j-1)*params.alpha;
        x_spoke = [x_hub x_hub + params.l*sin(spoke_angle)];
        y_spoke = [y_hub y_hub - params.l*cos(spoke_angle)];
        plot(x_spoke, y_spoke, 'b', 'LineWidth', 2);
        
        % Draw spoke endpoints
        plot(x_spoke(2), y_spoke(2), 'k.', 'MarkerSize', 10);
    end
    
    % Draw hub
    plot(x_hub, y_hub, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    % Draw contact point
    plot(x_contact, y_contact, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    
    axis equal;
    xlim([x_hub-2 x_hub+2]);

    title('Rimless Wheel Simulation');
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    drawnow;
   
end


%%

% Phase portrait plot
figure;
plot(Y(:,1), Y(:,2), 'b.', 'MarkerSize', 2);
hold on;

theta_range = linspace(min(Y(:,1)), max(Y(:,1)), 20);
omega_range = linspace(min(Y(:,2)), max(Y(:,2)), 20);
[THETA, OMEGA] = meshgrid(theta_range, omega_range);


dTHETA = OMEGA;
dOMEGA = params.g/params.l * sin(THETA);


mag = sqrt(dTHETA.^2 + dOMEGA.^2);
dTHETA = dTHETA./mag;
dOMEGA = dOMEGA./mag;


quiver(THETA, OMEGA, dTHETA, dOMEGA, 0.5, 'k', 'LineWidth', 0.5);


theta_line = linspace(min(Y(:,1)), max(Y(:,1)), 100);

theta_forward = -params.alpha/2 + params.gamma;
plot([theta_forward theta_forward], [min(Y(:,2)) max(Y(:,2))], 'g--', 'LineWidth', 1.5);

theta_backward = params.alpha/2 + params.gamma;
plot([theta_backward theta_backward], [min(Y(:,2)) max(Y(:,2))], 'r--', 'LineWidth', 1.5);

title('Phase Portrait (θ vs ω)');
xlabel('θ (rad)');
ylabel('ω (rad/s)');
grid on;


%%

% Dynamics function
function dydt = dynamics(t, y, params)
    theta = y(1);
    omega = y(2);
    
    % Equations of motion
    dydt = zeros(3,1);
    dydt(1) = omega;  % angular velocity
    dydt(2) = params.g/params.l * sin(theta);  % angular acceleration
    dydt(3) = 0;  % contact point doesn't move during continuous motion
end


function [value, isterminal, direction] = collisionEvent(t, y, params)
    theta = y(1);
    omega = y(2);  % angular velocity
    x_contact = y(3);
    
    % Calculate hub position
    y_contact = x_contact * tan(params.gamma);
    x_hub = x_contact - params.l*sin(theta);
    y_hub = y_contact + params.l*cos(theta);
    
    % Calculate height of all spoke endpoints relative to ground
    heights = zeros(params.n, 1);
    for i = 1:params.n
        spoke_angle = theta + (i-1)*params.alpha;
        x_end = x_hub + params.l*sin(spoke_angle);
        y_end = y_hub - params.l*cos(spoke_angle);
        y_ground = x_end * tan(params.gamma);
        heights(i) = y_end - y_ground;
    end
    
    % Find minimum height that's not the current stance leg
    [~, stance_leg] = min(abs(heights));
    if omega<0
        [~, stance_leg] = min((heights));
    end
    heights(stance_leg) = inf;  % exclude current stance leg
    
    % Find the minimum of remaining heights
    value = min(heights);
    
    
    
    direction=0;
    
    isterminal = 1;  % stop integration
end


function y_plus = collisionMap(y_minus, params)
    theta_minus = y_minus(1);
    omega_minus = y_minus(2);
    x_minus = y_minus(3);
    
    if omega_minus >= 0  % Forward motion (downhill)
        % Update contact point position forward
        x_plus = x_minus + 2*params.l*cos(params.gamma)*sin(params.alpha/2);
        
        % Reset angle to start of stance phase
        theta_plus = -params.alpha/2 + params.gamma;
        
        % Apply collision equation: ω⁺ = ω⁻cos(2π/n)
        omega_plus = omega_minus * cos(2*pi/params.n);
        
    else  % Backward motion (uphill)
        % Update contact point position backward
        x_plus = x_minus - 2*params.l*cos(params.gamma)*sin(params.alpha/2);
        
        % Reset angle for backward motion
        theta_plus = params.alpha/2 + params.gamma;
        
        % Apply collision equation with reversed direction
        omega_plus = omega_minus * cos(2*pi/params.n);
    end
    
    y_plus = [theta_plus; omega_plus; x_plus];
end

