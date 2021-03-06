% M-script for numerical integration of the attitude dynamics of a rigid 
% body represented by unit quaternions. The MSS m-files must be on your
% Matlab path in order to run the script.
%
% System:                      .
%                              q = T(q)w
%                              .
%                            I w - S(Iw)w = tau
% Control law:
%                            tau = -Kdw -kp(epsilon)
% 
% Definitions:             
%                            I = inertia matrix (3x3)
%                            S(w) = skew-symmetric matrix (3x3)
%                            T(q) = transformation matrix (4x3)
%                            tau = control input (3x1)
%                            w = angular velocity vector (3x1)
%                            q = unit quaternion vector (4x1)
%
% Author:                   2017-09-01 Maren Eidal & Helene H. Fossum

%% USER INPUTS
h = 0.1;                    % sample time (s)
N  = 2000;                  % number of samples

% model parameters
m = 100;                    % mass of the satellite [kg]
r = 2;                      % radius of the satellite [m]
I = m*r^2*diag([1 1 1]);    % inertia matrix
I_inv = inv(I);

%PD controller gains
kd = 20;
Kd = kd*eye(3);

kp = 1;

% constants
deg2rad = pi/180;   
rad2deg = 180/pi;

% initial Euler angles
phi = 10*deg2rad;           % roll
theta = -5*deg2rad;         % pitch     
psi = 15*deg2rad;           % yaw

q = euler2q(phi,theta,psi);   % transform initial Euler angles to q

w = [0 0 0]';                 % initial angular rates

table = zeros(N+1,14);        % memory allocation

%% FOR-END LOOP
for i = 1:N+1,
   t = (i-1)*h;                     % time
   tau  = -Kd*eye(3)*w - kp*q(2:4); % control law

   [phi,theta,psi] = q2euler(q);    % transform q to Euler angles
   [J,J1,J2] = quatern(q);          % kinematic transformation matrices
   
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
% plot in Euler angles
clf
figure(gcf)
subplot(311),plot(t,phi),xlabel('time (s)'),ylabel('deg'),title('Roll angle \phi'),grid
subplot(312),plot(t,theta),xlabel('time (s)'),ylabel('deg'),title('Pitch angle \theta'),grid
subplot(313),plot(t,psi),xlabel('time (s)'),ylabel('deg'),title('Yaw angle \psi'),grid
hold on 

figure
subplot(211),plot(t,w),xlabel('time (s)'),ylabel('deg/s'),title('Angular velocity w'),grid,
legend('p','q','r');
subplot(212),plot(t,tau),xlabel('time (s)'),ylabel('Nm'),title('Control input \tau'),grid,
legend('tau_1', 'tau_2', 'tau_3');
hold on 


%plot quaternion 
figure
subplot (211), plot(t,q(:,1)),xlabel('time(s)'),ylabel(''),title('\eta'),grid
subplot (212), plot(t,q(:,2)),xlabel('time(s)'),ylabel(''),title('\epsilon_1'),grid
hold on

figure
subplot (211), plot(t,q(:,3)),xlabel('time(s)'),ylabel(''),title('\epsilon_2'),grid
subplot (212), plot(t,q(:,4)),xlabel('time(s)'),ylabel(''),title('\epsilon_3'),grid
hold on 
