%% Synthetic INS simulation

% This framework was created by Barak Or.
% Version 0.1, released at July 6 2020.
% The Synthetic INS simulation provide the exact navigation solution based
% on the strapdown approch, with perfect IMU (no error for gyro and
% accelerometes).

% This code is dived into two main parts:
% 1) creating trajectory and solve the navigation equation.
% 2) Plotting the variuos graphes with ref. data structure, that can be
% implemented further in INS/DVL (or other fusion) environment.

%% General Setting
clc
close all
clear

set(0,'DefaultFigureWindowStyle','docked')
set(0,'defaulttextInterpreter','latex')

plotall=1;
scenario=1;
export_ref=0;
noise=0;
D2R = pi/180;     % deg.  to rad.
R2D = 180/pi;     % rad. to deg.

% Create time frame
dt    = 0.01; % sec.                        % IMU sampling interval
time  = 240;                                % simulation time.
ref.t = (0.01:dt:time)';                       % IMU time vector

%% Initialization
% ref.mat contains the reference data structure from which inertial
% sensors and DVL wil be simulated. It must contain the following fields:

%         t: Nx1 time vector (seconds).
%       lat: Nx1 latitude (radians).
%       lon: Nx1 longitude (radians).
%         h: Nx1 altitude (m).
%       vel: Nx3 NED velocities (m/s).
%      roll: Nx1 roll angles (radians).
%     pitch: Nx1 pitch angles (radians).
%       yaw: Nx1 yaw angle vector (radians).
%     DCMnb: Nx9 Direct Cosine Matrix nav-to-body. Each row contains
%            the elements of one DCM matrix ordered by columns as
%            [a11 a21 a31 a12 a22 a32 a13 a23 a33].
%      freq: sampling frequency (Hz).

%% Generate acc. and gyro in body frame
% Generate accelerometer measurements in body frame by setting fb, and
% angular velocity measurements in body frame by setting wb. For simplicity,
% they are initialized with zero (expect the 3rd value in fb).

%% INS only

% Constant matrices
I = eye(3);
O = zeros(3);
count=1;
L = length(ref.t);

% Attitude
roll_e  = zeros (L, 1);
pitch_e = zeros (L, 1);
yaw_e   = zeros (L, 1);

DCMnb = euler2dcm([roll_e(1); pitch_e(1); yaw_e(1);]);
DCMbn = DCMnb';
qua   = euler2qua([roll_e(1) pitch_e(1) yaw_e(1)]);
DCMnb_store(1,:)=[1 0 0 0 1 0 0 0 1];
% Initialize velocity vector at INS (N.E.D)
vel_e(1,:) = [5 0 0];
wb =[0 0 0];
wb_store(1,:) =wb;
fb_store(1,:)=[0 0 9.8];
v_total(1)=norm([vel_e],2);

% Initialize position vector at INS
h_e(1)   = -5; %m
lat_e(1) = pi/3; %rad
lon_e(1) = pi/3; %rad
t=0.01;
t_store(1)=0;


for i=2:L
    
    if scenario == 1  % circular diving
        
        fb_new(i,:)=[0 0 0]';
        wb(i,:) =[0 0 0];
        
        if t>50 && t< 60
            fb_new(i,:) = [-0.5 0.5 0]';
            
        end
        
        if t>110 && t< 120
            fb_new(i,:) = [ -0.5 -0.5 0]';
        end
        
        if t>170 && t< 180
            fb_new(i,:) = [ 0.5 -0.5 0]';
        end
        
    end
    
    % Turn-rates update
    omega_ie_n = earthrate(lat_e(i-1));
    omega_en_n = transportrate(lat_e(i-1), vel_e(i-1,1), vel_e(i-1,2), h_e(i-1));
    
    % In order to obtain perfect IMU readings, the angular velocity is
    % exactly as shuld be measured in body frame.
    wb_store(i,:) = wb(i,:);
    
    % Attitude update
    [qua_n, DCMbn, euler] = att_update(wb(i,:), DCMbn, qua, ...
        omega_ie_n, omega_en_n, dt, 'dcm');
    roll_e(i) = euler(1);
    pitch_e(i)= euler(2);
    yaw_e(i)  = euler(3);
    qua = qua_n;
    DCMnb= DCMbn';
    DCMnbstore(i,:)=DCMnb(:);
    
    
    % Gravity update
    gn = gravity(lat_e(i-1), h_e(i-1));
    fb_clean = DCMnb * gn';
    % In order to obtain perfect IMU readings, the gravity vector is
    % substituted in place of fb.
    fb=fb_new(i,:)'+fb_clean +0*0.1*randn(3,1);
    fb_store(i,:) = fb';
    fn = (DCMbn * fb);
    
    % Velocity update
    vel_n = vel_update(fn, vel_e(i-1,:), omega_ie_n, omega_en_n, gn', dt);
    vel_e (i,:) = vel_n;
    v_store(i,:)=vel_e (i,:);
    v_total(i)=norm([vel_n],2);
    
    % Position update
    pos = pos_update([lat_e(i-1) lon_e(i-1) h_e(i-1)], vel_e(i,:), dt);
    
    lat_e(i)   = pos(1);
    lon_e(i)   = pos(2);
    h_e(i)     = pos(3);
    
    t_store(i) = t_store(i-1)+dt;
    t=t+dt;
end

%% Creating reference (ground truth) for further analysis.
wb_store(1,:) =wb_store(2,:);

ref.lat   = lat_e';
ref.lon   = lon_e';
ref.h     = h_e';
ref.v     = v_store;
ref.roll  = roll_e;
ref.pitch = pitch_e;
ref.yaw   = yaw_e;
ref.DCMnb = DCMnbstore;
ref.fb    = fb_store;
ref.wb    = wb_store;
ref.freq=100;
ref.t = t_store;
if export_ref==1
    save ref_clean.mat;
end
%% Plot
if plotall==1
    lat_e=lat_e.*R2D;
    lon_e=lon_e.*R2D;
    
    figure
    set(gca,'Fontsize',12);
    plot3(lat_e,lon_e,h_e)
    zlim([-20,-2])
    grid minor
    xlabel('Longtitude $\lambda$[deg]')
    ylabel('Latitude $\phi$ [deg]')
    zlabel('Height, $h$ [m]')
    
    figure
    subplot(2,1,1)
    plot(lat_e,lon_e)
    title('Latitude,  $\phi$ vs. Longtitude, $\lambda$')
    xlabel('Longtitude $\lambda$[deg]')
    ylabel('Latitude $\phi$ [deg]')
    grid minor
    subplot(2,1,2)
    plot(t_store,h_e)
    xlim([0,time])
    ylim([-6,-4])
    if scenario == 1
        ylim([-10,-4])
    end
    title('Height, $h$ vs. Time')
    xlabel('Time [s]')
    ylabel('h [m]')
    grid minor
    
    figure
    subplot(2,1,1)
    plot(t_store,v_store (:,1),t_store,v_store (:,2),t_store,v_store (:,3))
    title('Velocity')
    xlim([0,time])
    xlabel('Time [s]')
    legend('$v^N$','$v^E$','$v^D$','Interpreter','latex')
    grid minor
    
    subplot(2,1,2)
    plot(t_store,v_total)
    xlim([0,time])
    ylim([3,6])
    xlabel('Time [s]')
    grid minor
    title('Speed ($\left\| {{v_j}} \right\|_2^2)$')
    
    figure
    subplot(3,1,1)
    plot(t_store,roll_e)
    ylabel('$\phi [rad]$')
    xlim([0,time])
    xlabel('Time [s]')
    grid minor
    subplot(3,1,2)
    plot(t_store,pitch_e)
    ylabel('$\theta$')
    xlim([0,time])
    xlabel('Time [s]')
    grid minor
    subplot(3,1,3)
    plot(t_store,yaw_e)
    xlim([0,time])
    ylabel('$\psi$')
    xlabel('Time [s]')
    grid minor
    
    
    figure
    title('IMU Readings')
    subplot(2,1,1)
    plot(t_store,wb_store(:,1),t_store,wb_store(:,2),t_store,wb_store(:,3))
    xlim([0,time])
    legend('$\omega_b^x$','$\omega_b^y$','$\omega_b^z$','Interpreter','latex')
    grid minor
    xlabel('Time [s]')
    ylabel('$\omega_b^i [rad/s]$')
    subplot(2,1,2)
    plot(t_store,fb_store(:,1),t_store,fb_store(:,2),t_store,fb_store(:,3))
    ylim([-2,12])
    xlim([0,time])
    legend('$f_b^x$','$f_b^y$','$f_b^z$','Interpreter','latex')
    grid minor
    xlabel('Time [s]')
    ylabel('$f_b^i [m/s^2]$')
end