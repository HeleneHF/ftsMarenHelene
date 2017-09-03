% M-script for numerical integration of the attitude dynamics of a rigid 
% body represented by unit quaternions. The MSS m-files must be on your
% Matlab path in order to run the script.
%
% System:                      .
%                              q = T(q)w
%                              .
%                            I w - S(Iw)w = tau
% Control law:
%                            tau = -Kd(w_tilde) -kp(epsilon_tilde)
% 
% Definitions:             
%                            I = inertia matrix (3x3)
%                            S(w) = skew-symmetric matrix (3x3)
%                            T(q) = transformation matrix (4x3)
%                            tau = control input (3x1)
%                            w = angular velocity vector (3x1)
%                            q = unit quaternion vector (4x1)
%
% Author:                   2017-09-03 Maren Eidal & Helene H. Fossum


%% USER INPUTS
h = 0.1;                    % sample time (s)
N  = 8000;                  % number of samples

% model parameters
m = 100;                    % mass of the satellite [kg]
r = 2;                      % radius of the satellite [m]
I = m*r^2*diag([1 1 1]);    % inertia matrix
I_inv = inv(I);

%PD controller gains
kd = 300;
Kd = kd*eye(3);

kp = 10;

% constants
deg2rad = pi/180;   
rad2deg = 180/pi;

phi = 10*deg2rad;           % initial Euler angles
theta = -5*deg2rad;
psi = 15*deg2rad;

q = euler2q(phi,theta,psi);   % transform initial Euler angles to q

w = [0 0 0]';                 % initial angular rates
w_d = [0 0 0]';
table = zeros(N+1,24);        % memory allocation

%% FOR-END LOOP
for i = 1:N+1
   t = (i-1)*h;                     % time
   
    % Actual position
    eta = q(1);
    e = q(2:4);
    
    % Desired position in Euler angles
    phi_d = (10*sin(0.1*t))*deg2rad;
    theta_d = 0;
    psi_d = (15*cos(0.05*t))*deg2rad; 

    %Til Maren
%     phi_d = 0;
%     psi_d = 0;
   
    % Desired position in unit quaternions
    q_d = euler2q(phi_d,theta_d,psi_d);
    eta_d = q_d(1);
    e_d = q_d(2:4);

    % w_tilde
    T_inv = [1 0 -sin(theta_d); 0 cos(phi_d) cos(theta_d)*sin(phi_d); 0 -sin(phi_d) cos(theta_d)*cos(phi_d)]; 
    Theta_big = [phi_d,theta_d,psi_d];
    Theta_big_dot = [deg2rad*cos(0.1*t); 0; deg2rad*(-0.75)*sin(0.05*t)];
    w_d = T_inv*Theta_big_dot;

    % The quaternion error using the conjugate of a quaternion
    q_tilde = [eta_d*eta + e_d'*e; eta_d*e - eta*e_d + Smtrx(-e_d)*e];
    eta_tilde = q_tilde(1);
    e_tilde = q_tilde(2:4); 
    
    w_tilde = w - w_d;
    
   tau  = -Kd*eye(3)*w_tilde - kp*e_tilde; % control law

   [phi,theta,psi] = q2euler(q);    % transform q to Euler angles
   [J,J1,J2] = quatern(q);          % kinematic transformation matrices
   
   % Error from desired to actual pos. in Euler angles
   phi_error = phi - phi_d;
   theta_error = theta - theta_d;
   psi_error = psi - psi_d; 
   
   q_dot = J2*w;                        % quaternion kinematics
   w_dot = I_inv*(Smtrx(I*w)*w + tau);  % rigid-body kinetics
   
   table(i,:) = [t q' phi theta psi w' tau' phi_error, theta_error, psi_error, phi_d,theta_d,psi_d,eta_tilde,w_d'];  % store data in table
   
   q = q + h*q_dot;	             % Euler integration
   w = w + h*w_dot;
   
   q  = q/norm(q);               % unit quaternion normalization
end 

%% Storing data
t       = table(:,1);  
q       = table(:,2:5); 

phi     = rad2deg*table(:,6);
theta   = rad2deg*table(:,7);
psi     = rad2deg*table(:,8);

w       = rad2deg*table(:,9:11);  
tau     = table(:,12:14);

phi_error = rad2deg*table(:,15);
theta_error = rad2deg*table(:,16);
psi_error = rad2deg*table(:,17);

phi_d = rad2deg*table(:,18);
theta_d = rad2deg*table(:,19);
psi_d = rad2deg*table(:,20);

phi_a = [phi,phi_d];
theta_a = [theta, theta_d];
psi_a = [psi, psi_d];

eta_tilde = table(:,21);
w_d = rad2deg*table(:,22:24);

%% Plot
% plot in Euler angles
close all

clf
figure(1) 
subplot(311),plot(t,phi_a),xlabel('time (s)'),ylabel('deg'),title('Roll angle \phi'),grid, legend('phi','phi_d')
subplot(312),plot(t,theta_a),xlabel('time (s)'),ylabel('deg'),title('Pitch angle \theta'),grid, legend('theta','theta_d')
subplot(313),plot(t,psi_a),xlabel('time (s)'),ylabel('deg'),title('Yaw angle \psi'),grid, legend('psi','psi_d')
hold on 

figure(2)
subplot(211),plot(t,w),xlabel('time (s)'),ylabel('deg/s'),title('Angular velocity w'),grid,
legend('p','q','r');
subplot(212),plot(t,tau),xlabel('time (s)'),ylabel('Nm'),title('Control input \tau'),grid,
legend('tau_1', 'tau_2', 'tau_3');
hold on 

 %plot quaternion 
%  figure(3)
%  subplot (211), plot(t,q(:,1)),xlabel('time(s)'),ylabel(''),title('\eta'),grid
%  subplot (212), plot(t,q(:,2)),xlabel('time(s)'),ylabel(''),title('\epsilon_1'),grid
%  hold on
% 
% figure(4)
% subplot (211), plot(t,q(:,3)),xlabel('time(s)'),ylabel(''),title('\epsilon_2'),grid
% subplot (212), plot(t,q(:,4)),xlabel('time(s)'),ylabel(''),title('\epsilon_3'),grid
% hold on 

% Plot error in euler angles
figure(5)
subplot(311),plot(t,phi_error),xlabel('time (s)'),ylabel('deg'),title('Error roll angle \phi'),grid
subplot(312),plot(t,theta_error),xlabel('time (s)'),ylabel('deg'),title('Error pitch angle \theta'),grid
subplot(313),plot(t,psi_error),xlabel('time (s)'),ylabel('deg'),title('Error yaw angle \psi'),grid
hold on 

% figure(6)
% plot(t,eta_tilde),xlabel('time(s)'),ylabel(''),title('\eta~'),grid

% figure(7)
% plot(t,w_d);