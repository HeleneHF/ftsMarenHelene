%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name:         Attitude 4                                          %
% Assignment:   Assignment 1 Part 2                                 %
% Author:       2017-09-13 Helene H.Fossum                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% USER INPUTS
% simulation parameters
h = 0.1;                    % sample time (s)
N  = 8000;                  % number of samples
current = false;             % turns the current on/off

% constants
deg2rad = pi/180;   
rad2deg = 180/pi;

% model parameters
p = [0;0;10];   % initial position, 10 [m] depth
Radius = 100;   % [m] radius of the circle
U = 1.5;        % [m/s] speed without current
U_c = 0.6;      % [m/s] speed of the ocean current
omega = 0.015;  % from problem 2.2

phi = 0*deg2rad;    % [rad] roll angle
theta = 2*deg2rad;  % [rad] pitch angle
psi = 30*deg2rad;   % [rad] yaw angle

alpha_c = 10*deg2rad; % [rad]
beta_c = 45*deg2rad;  % [rad]  

% rotation matrix FLOW to NED
Rot_flow_1 = [cos(beta_c)*cos(alpha_c), sin(beta_c), cos(beta_c)*sin(alpha_c)];
Rot_flow_2 = [-sin(beta_c)*cos(alpha_c), cos(beta_c), -sin(beta_c)*sin(alpha_c)];
Rot_flow_3 = [-sin(alpha_c), 0, cos(alpha_c)];

Rot_flow = [Rot_flow_1; Rot_flow_2; Rot_flow_3]; % final

% rotation matrix NED to BODY
[J,J11,J22] = eulerang(phi,theta,psi); 

table = zeros(N+1,11);        % memory allocation

%% FOR-END LOOP
for i = 1:N+1   
   t = (i-1)*h; % time
    
    % velocites circle no current
    u = U*cos(omega*t); % [m/s] undirsturbed velocity x-dir (body)
    v = U*sin(omega*t); % [m/s] undirsturbed velocity y-dir (body)
    w = 0;              % [m/s] undirsturbed velocity z-dir (body)
    
    
    if (current) % current on/off
    % ocean current velocities
        u_c = U_c*cos(omega*t);    % [m/s] ocean current velocity x-dir (flow)
        v_c = U_c*sin(omega*t);    % [m/s] ocean current velocity y-dir (flow)   
        w_c = 0;                   % [m/s] ocean current velocity z-dir (flow)
    else 
        u_c = 0;    % [m/s] ocean current velocity x-dir (flow)
        v_c = 0;    % [m/s] ocean current velocity y-dir (flow)   
        w_c = 0;    % [m/s] ocean current velocity z-dir (flow)
    end
    
    % for problem 2.3
    v_curr_flow = [u_c; v_c; w_c];              % [m/s] ocean current velocities collected (flow)
    v_curr_body = Rot_flow'*J11*v_curr_flow;    % [m/s] ocean current velocities collected (body)
    
    % relative velocities
    u_r = u - u_c;    % [m/s] relative velocity x-dir (flow)
    v_r = v - v_c;    % [m/s] relative velocity y-dir (flow)   
    w_r = w - w_c;    % [m/s] realtive velocity z-dir (flow)    
    
    v_rel_flow = [u_r; v_r; w_r]; % [m/s] relative velocities collected (flow)
    v_rel_body = Rot_flow'*J11*v_rel; % [m/s] relative velocities collected (body)
    
    
    % NB!! NO PROPER EXPRESSIONS
    beta_r = asin(v_r/u_r); % [rad] sideslip angle     
    beta = asin(v/U);       % [rad] crab angle
    chi = beta + psi;       % [rad] course angle
    
    U_r = sqrt(u_r^2 + v_r^2 + w_r^2); % [m/s] velocity
    
    % equations of motion 
    p_dot = J11*v_rel_flow; % eq. 2.40 EoM
    
    table(i,:) = [t,u_r,v_r,w_r,beta_r,beta,chi,U_r,p'];  % store data in table
    
    p = p + h*p_dot; % euler integration
   
end 

%% Defining variables
t = table(:,1);

u_r = table(:,2);    % [m/s] relative velocity x-dir (body?)
v_r = table(:,3);    % [m/s] relative velocity y-dir (body?)   
w_r = table(:,4);    % [m/s] realtive velocity z-dir (body?)       
    
beta_r = table(:,5)*rad2deg; % [degree] sideslip angle     
beta = table(:,6)*rad2deg;   % [degree] crab angle
chi = table(:,7)*rad2deg;    % [degree] course angle

U_r = table(:,8); % [m/s] velocity

N = table(:,9); % position north
E = table(:,10); % position east   
D = table(:,11); % position down   

%% Plot
close all

clf
figure(1) % postion in NED-frame
plot(E,N),xlabel('East'),ylabel('North'), title ('Postion in NED'), grid 
hold on

figure(2) % velocity components
subplot(311),plot(t,u_r),xlabel('time[s]'),ylabel('u_r'), title ('Relative velocity x-dir'), grid
subplot(312),plot(t,v_r),xlabel('time[s]'),ylabel('v_r'), title ('Relative velocity y-dir'), grid
subplot(313),plot(t,w_r),xlabel('time[s]'),ylabel('w_r'), title ('Relative velocity z-dir'), grid
hold on 


figure(3) % velocity
if (current)
    plot(t,U_r),xlabel('time[s]'),ylabel('velocity [m/s]'), title('Relative velocity U_r'), grid
else
   plot(t,U_r),xlabel('time[s]'),ylabel('velocity [m/s]'), title('Undisturbed velocity U'), grid 
end
hold on 

figure(4) % sideslip, crab and course angle 
subplot(311),plot(t,beta_r),xlabel('time[s]'),ylabel('angle[deg]'), title ('Sideslip angle \beta_r'), grid
subplot(312),plot(t,beta),xlabel('time[s]'),ylabel('angle[deg]'), title ('Crab angle \beta'), grid
subplot(313),plot(t,chi),xlabel('time[s]'),ylabel('angle[deg]'), title ('Course angle \chi'), grid
hold on 





