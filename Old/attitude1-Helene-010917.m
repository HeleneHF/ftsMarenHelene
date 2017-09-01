



% LARGE MISTAKE IN I MATRIX, DO NOT USE



% M-script for numerical integration of the attitude dynamics of a rigid 
% body represented by unit quaternions. The MSS m-files must be on your
% Matlab path in order to run the script.
%
% System:                      .
%                              q = T(q)w
%                              .
%                            I w - S(Iw)w = tau
% Control law:
%                            tau = -Kdw -kp(epsilon_tau)
% 
% Definitions:             
%                            I = inertia matrix (3x3)
%                            S(w) = skew-symmetric matrix (3x3)
%                            T(q) = transformation matrix (4x3)
%                            tau = control input (3x1)
%                            w = angular velocity vector (3x1)
%                            q = unit quaternion vector (4x1)
%
% Author:                   2017-08-28 Maren Eidal Helene H. Fossum  

%% USER INPUTS
h = 0.1;                     % sample time (s)
N  = 5000;                   % number of samples

% model parameters
m = 100;                    % mass of the satellite [kg]
r = 2;                      % radius of the satellite [m]
I = diag([m*r^2 m*r^2 m*r^2]);    % inertia matrix (sphere)
I_inv = inv(I);

% PD controller gains
kp = 1;     Kp = kp*I(3);
kd = 20;    Kd = kd*I(3);

% constants
deg2rad = pi/180;   
rad2deg = 180/pi;

% Initial angles (euler)
phi = -10*deg2rad;            % initial Euler angles
theta = 10*deg2rad;
psi = 5*deg2rad;

q = euler2q(phi,theta,psi);   % transform initial Euler angles to q

w = [0 0 0]';                 % initial angular rates



table = zeros(N+1,14);        % memory allocation

%% FOR-END LOOP
for i = 1:N+1,
   t = (i-1)*h;                  % time
   
   e = q(2:4);                   % Imaginary parts of the u. quaternion
   
   tau = -Kd*w -kp*e;            % control law (Prob. 1.2)

   [phi,theta,psi] = q2euler(q); % transform q to Euler angles
   [J,J1,J2] = quatern(q);       % kinematic transformation matrices
   
   q_dot = J2*w;                        % quaternion kinematics
   w_dot = I_inv*(Smtrx(I*w)*w + tau);  % rigid-body kinetics
   
   table(i,:) = [t q' phi theta psi w' tau'];  % store data in table
   
   q = q + h*q_dot;	             % Euler integration
   w = w + h*w_dot;
   
   q  = q/norm(q);               % unit quaternion normalization
end 

%% PLOT FIGURES
t       = table(:,1);  
q       = table(:,2:5); 
phi     = rad2deg*table(:,6);
theta   = rad2deg*table(:,7);
psi     = rad2deg*table(:,8);
w       = rad2deg*table(:,9:11);  
tau     = table(:,12:14);

%% Plot
% % clf
 
figure(gcf)
subplot(211);plot(t,phi),xlabel('time (s)'),ylabel('deg'),title('Roll angle \phi'),grid
subplot(212);plot(t,theta),xlabel('time (s)'),ylabel('deg'),title('Pitch angle \theta'),grid
hold on 
 
figure
subplot(311),plot(t,psi),xlabel('time (s)'),ylabel('deg'),title('Yaw angle \psi'),grid
subplot(312),plot(t,w),xlabel('time (s)'),ylabel('deg/s'),title('Angular velocity w'),grid
subplot(313),plot(t,tau),xlabel('time (s)'),ylabel('Nm'),title('Control input \tau'),grid
