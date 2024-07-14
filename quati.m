clc
clear
close all

% MATLAB script to visualize quaternion components over time

%% Define initial conditions %%
s0 = [0; 0; 0]; % initial position (m)
v0 = [-1; 0; 2]; % initial velocity (m/s)
omega0 = [0.5; -0.7; 0.3]; % initial angular velocity (rad/s)
q_bar0 = [1; 0; 0; 0]; % initial quaternion (unit quaternion)

% Define inertia matrix (assuming a cuboid for simplicity)
mass = 1; % kg
a = 0.1; b = 0.1; c = 0.1; % dimensions of the cuboid (m)
Ixx = mass/12 * (b^2 + c^2);
Iyy = mass/12 * (a^2 + c^2);
Izz = mass/12 * (a^2 + b^2);
I_bar = diag([Ixx, Iyy, Izz]);

% Define simulation parameters
%tspan = [0 20]; % time span for the simulation (s)
tspan = linspace(0,20,500);
y0 = [s0; v0; q_bar0; omega0]; % initial state vector

%% Solve the differential equations %%
tol = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, y] = ode45(@(t, y) rigid_body_dynamics(t, y, I_bar), tspan, y0, tol);

% Extract the quaternion components over time
q_bar = y(:, 7:10);

%% Plot quaternion components vs. time %%
figure;
subplot(4, 1, 1);
plot(t, q_bar(:, 1), 'r', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('q0');
title('Quaternion Components vs. Time');

subplot(4, 1, 2);
plot(t, q_bar(:, 2), 'g', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('q1');

subplot(4, 1, 3);
plot(t, q_bar(:, 3), 'b', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('q2');

subplot(4, 1, 4);
plot(t, q_bar(:, 4), 'k', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('q3');

eulerAngles = quaternionsToEulerAngles(q_bar, t);

plotEulerAngleBodyCone(eulerAngles, t)

animate_cube(eulerAngles, t);

% Define the differential equations
function dydt = rigid_body_dynamics(t, y, I_bar)
    % Extract variables from state vector
    s = y(1:3);
    v = y(4:6);
    q_bar = y(7:10);
    omega = y(11:13);

    % Translational motion (assuming no external forces)
    dsdt = v;
    dvdt = [0; 0; 0];

    % Rotational motion
    q = q_bar';
    omega_quat = [0; omega];
    dqdt = 0.5 * quatmultiply(q, omega_quat');
    dqdt = dqdt';

    domega = I_bar \ (-cross(omega, I_bar * omega));

    % Combine derivatives into a single vector
    dydt = [dsdt; dvdt; dqdt; domega];
end

function eulerAngles = quaternionsToEulerAngles(quaternion, t)
    % Convert quaternions to Euler angles (yaw, pitch, roll)
    % Inputs:
    %   quaternion - Nx4 matrix of quaternions
    %   t - Nx1 time vector
    % Outputs:
    %   eulerAngles - Nx3 matrix of Euler angles [yaw, pitch, roll]

    % Initialize output matrix
    eulerAngles = zeros(length(t), 3);

    for i = 1:length(t)
        q0 = quaternion(i, 1);
        q1 = quaternion(i, 2);
        q2 = quaternion(i, 3);
        q3 = quaternion(i, 4);

        % Calculate Euler angles
        yaw = atan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2^2 + q3^2));
        pitch = asin(2 * (q0 * q2 - q3 * q1));
        roll = atan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1^2 + q2^2));

        eulerAngles(i, :) = [yaw, pitch, roll];
    end
end

function plotEulerAngleBodyCone2(eulerAngles, t)
    % Function to plot the 3D Euler angle body cone
    % Inputs:
    %   eulerAngles - Nx3 matrix of Euler angles [yaw, pitch, roll]
    %   t - Nx1 time vector

    % Define the length of the cone lines
    cone_length = 1;

    % Initialize arrays to store the tip positions of the cone lines
    X = zeros(length(t), 1);
    Y = zeros(length(t), 1);
    Z = zeros(length(t), 1);

    for i = 1:length(t)
        yaw = eulerAngles(i, 1);
        pitch = eulerAngles(i, 2);
        roll = eulerAngles(i, 3);

        % Compute the rotation matrix from the Euler angles
        Rz = [cos(yaw) -sin(yaw) 0; sin(yaw) cos(yaw) 0; 0 0 1]; % Yaw
        Ry = [cos(pitch) 0 sin(pitch); 0 1 0; -sin(pitch) 0 cos(pitch)]; % Pitch
        Rx = [1 0 0; 0 cos(roll) -sin(roll); 0 sin(roll) cos(roll)]; % Roll
        R = Rz * Ry * Rx;

        % Define the body-fixed axis (e.g., the z-axis of the body frame)
        body_axis = [0; 0; cone_length];

        % Rotate the body axis to get its orientation in the inertial frame
        rotated_axis = R * body_axis;

        % Store the tip position of the cone line
        X(i) = rotated_axis(1);
        Y(i) = rotated_axis(2);
        Z(i) = rotated_axis(3);
    end

    % Plot the body cone
    figure;
    plot3(X, Y, Z, 'LineWidth', 2);
    hold on;
    plot3([0 X(1)], [0 Y(1)], [0 Z(1)], 'r', 'LineWidth', 2); % Initial orientation
    plot3([0 X(end)], [0 Y(end)], [0 Z(end)], 'b', 'LineWidth', 2); % Final orientation
    grid on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('3D Euler Angle Body Cone');
    axis equal;
end

function plotEulerAngleBodyCone(eulerAngles, t)
    % Function to plot the 3D Euler angle body cone and space cone
    % Inputs:
    %   eulerAngles - Nx3 matrix of Euler angles [yaw, pitch, roll]
    %   t - Nx1 time vector

    % Define the length of the cone lines
    cone_length = 1;

    % Initialize arrays to store the tip positions of the cone lines for body cone
    X_body = zeros(length(t), 1);
    Y_body = zeros(length(t), 1);
    Z_body = zeros(length(t), 1);
    
    % Initialize arrays to store the tip positions of the cone lines for space cone
    X_space = zeros(length(t), 1);
    Y_space = zeros(length(t), 1);
    Z_space = zeros(length(t), 1);

    % Initialize the inertial frame fixed axis for space cone
    inertial_axis = [0; 0; cone_length];

    for i = 1:length(t)
        yaw = eulerAngles(i, 1);
        pitch = eulerAngles(i, 2);
        roll = eulerAngles(i, 3);

        % Compute the rotation matrix from the Euler angles
        Rz = [cos(yaw) -sin(yaw) 0; sin(yaw) cos(yaw) 0; 0 0 1]; % Yaw
        Ry = [cos(pitch) 0 sin(pitch); 0 1 0; -sin(pitch) 0 cos(pitch)]; % Pitch
        Rx = [1 0 0; 0 cos(roll) -sin(roll); 0 sin(roll) cos(roll)]; % Roll
        R = Rz * Ry * Rx;

        % Define the body-fixed axis (e.g., the z-axis of the body frame)
        body_axis = [0; 0; cone_length];

        % Rotate the body axis to get its orientation in the inertial frame
        rotated_body_axis = R * body_axis;

        % Rotate the inertial axis to get its orientation in the body frame
        rotated_inertial_axis = R' * inertial_axis;

        % Store the tip position of the cone lines
        X_body(i) = rotated_body_axis(1);
        Y_body(i) = rotated_body_axis(2);
        Z_body(i) = rotated_body_axis(3);

        X_space(i) = rotated_inertial_axis(1);
        Y_space(i) = rotated_inertial_axis(2);
        Z_space(i) = rotated_inertial_axis(3);
    end

    % Plot the body cone
    figure;
    plot3(X_body, Y_body, Z_body, 'LineWidth', 2);
    hold on;
    plot3([0 X_body(1)], [0 Y_body(1)], [0 Z_body(1)], 'r', 'LineWidth', 2); % Initial orientation
    plot3([0 X_body(end)], [0 Y_body(end)], [0 Z_body(end)], 'b', 'LineWidth', 2); % Final orientation
    
    % Plot the space cone
    plot3(X_space, Y_space, Z_space, 'g--', 'LineWidth', 2);
    plot3([0 X_space(1)], [0 Y_space(1)], [0 Z_space(1)], 'm--', 'LineWidth', 2); % Initial orientation
    plot3([0 X_space(end)], [0 Y_space(end)], [0 Z_space(end)], 'c--', 'LineWidth', 2); % Final orientation

    % Formatting the plot
    grid on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('3D Euler Angle Body Cone and Space Cone');
    axis equal;
    legend('Body Cone', 'Initial Body Orientation', 'Final Body Orientation', 'Space Cone', 'Initial Space Orientation', 'Final Space Orientation');
    hold off;
end

function animate_cube(euler_angles, time_intervals)
    % Function to animate a 3D cube rotating based on given Euler angles and time intervals

    % Create a figure and axis for the animation
    figure;
    ax = axes('XLim', [-2 2], 'YLim', [-2 2], 'ZLim', [-2 2]);
    view(3);
    grid on;
    hold on;

    % Define the vertices of a cube
    vertices = [-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1];

    % Define the faces of the cube using the vertices
    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];

    % Plot the initial cube
    cube = patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'cyan', 'FaceAlpha', 0.3, 'EdgeColor', 'r');
    axis square

    % Number of frames
    num_frames = size(euler_angles, 1);

    % Animation loop
    for k = 1:num_frames
        % Get the current Euler angles
        theta = euler_angles(k, :);

        % Compute the rotation matrix
        R = euler_to_rotation_matrix(theta);

        % Rotate the vertices
        rotated_vertices = (R * vertices')';

        % Update the cube vertices
        set(cube, 'Vertices', rotated_vertices);

        % Pause for the specified time interval
        pause(time_intervals(k)/8);
    end
end

function R = euler_to_rotation_matrix(theta)
    % Function to compute the rotation matrix from Euler angles
    % theta: [yaw, pitch, roll]

    yaw = theta(1);
    pitch = theta(2);
    roll = theta(3);

    R_yaw = [cos(yaw), -sin(yaw), 0; sin(yaw), cos(yaw), 0; 0, 0, 1];
    R_pitch = [cos(pitch), 0, sin(pitch); 0, 1, 0; -sin(pitch), 0, cos(pitch)];
    R_roll = [1, 0, 0; 0, cos(roll), -sin(roll); 0, sin(roll), cos(roll)];

    R = R_yaw * R_pitch * R_roll;
end

