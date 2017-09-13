%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name:         Attitude 4                                          %
% Assignment:   Assignment 1 Part 2                                 %
% Author:       2017-09-13 Helene H.Fossum                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% USER INPUTS
% simulation parameters
h = 0.1;                    % sample time (s)
N  = 4000;                  % number of samples

% constants
deg2rad = pi/180;   
rad2deg = 180/pi;

% model parameters
p_init = [0,0,0]; % starts in the origin
R = 100;    % [m] radius of the circle
U = 1.5;    % [m/s] speed without current
U_c = 0.6;  % [m/s] speed of the ocean current
omega = 0.015; % from problem 2.2

phi = 0;    % [degree] roll angle
tetha = 2;  % [degree] pitch angle
psi = 30;   % [degree] yaw angle

alpha_c = 10*deg2rad; % [rad]
beta_c = 45*deg2rad;  % [rad]  

table = zeros(N+1,10);        % memory allocation

%% FOR-END LOOP
for i = 1:N+1
   t = (i-1)*h;                     % time
    
   % position
    N = 1; % north position 
    E = 2; % east position
    
    % velocites circle no current
    u = U*cos(omega*t); % [m/s] undirsturbed velocity x-dir (body?)
    v = U*sin(omega*t); % [m/s] undirsturbed velocity y-dir (body?)
    w = 0;              % [m/s] undirsturbed velocity z-dir (body?)
 
    % ocean current velocities
    u_c = 1;    % [m/s] ocean current velocity x-dir (body?)
    v_c = 2;    % [m/s] ocean current velocity y-dir (body?)   
    w_c = 3;    % [m/s] ocean current velocity z-dir (body?)    
   
    % retive velocities
    u_r = u - u_c;    % [m/s] relative velocity x-dir (body?)
    v_r = v - v_c;    % [m/s] relative velocity y-dir (body?)   
    w_r = w - w_c;    % [m/s] realtive velocity z-dir (body?)    
    
    % NB!! NO PROPER EXPRESSIONS
    beta_r = asin(v_r/u_r); % [rad] sideslip angle     
    beta = asin(v/U);       % [rad] crab angle
    chi = beta_r + yaw*deg2rad; % [rad] course angle
    
    U_r = sqrt(u_r^2 + v_r^2 + w_r^2); % [m/s] velocity
    
    table(i,:) = [t,N,E,u_r,v_r,w_r,beta_r,beta,chi,U_r];  % store data in table
   
end 

%% Defining variables
t = table(:,1);

N = table(:,2); % north position 
E = table(:,3); % east position

u_r = table(:,4);    % [m/s] relative velocity x-dir (body?)
v_r = table(:,5);    % [m/s] relative velocity y-dir (body?)   
w_r = table(:,6);    % [m/s] realtive velocity z-dir (body?)       
    
beta_r = table(:,7)*rad2deg; % [degree] sideslip angle     
beta = table(:,8)*rad2deg;   % [degree] crab angle
chi = table(:,9)*rad2deg;    % [degree] course angle

U_r = table(:,10); % [m/s] velocity

%% Plot
close all

clf
figure(1) % postion in NED-frame
plot(N,E),xlabel('East'),ylabel('North'), title ('Postion in NED'), grid 
hold on

figure(2)
subplot(311),plot(t,u_r),xlabel('time[s]'),ylabel('u_r'), title ('Relative velocity x-dir'), grid
subplot(312),plot(t,v_r),xlabel('time[s]'),ylabel('v_r'), title ('Relative velocity y-dir'), grid
subplot(313),plot(t,w_r),xlabel('time[s]'),ylabel('w_r'), title ('Relative velocity z-dir'), grid
hold on 

figure(3)
plot(t,U_r),xlabel('time[s]'),ylabel('velocity [m/s]'), title('Relative velocity U_r'), grid
hold on

% NB! SOMETHING WRONG, CHI MISSING!
figure(4) % sideslip, crab and course angle 
plot(t,beta_r,beta,chi), xlabel('time [s]'),ylabel('deg'), title('Sideslip, crab and course angle'), grid, legend('\beta_r', '\beta', '\chi') 
hold on 





