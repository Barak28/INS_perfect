% load ref
load ref.mat
ref.freq=100;


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

%% ADIS16405 IMU error profile

% IMU data structure:
%         t: Ix1 time vector (seconds).
%        fb: Ix3 accelerations vector in body frame XYZ (m/s^2).
%        wb: Ix3 turn rates vector in body frame XYZ (radians/s).
%       arw: 1x3 angle random walks (rad/s/root-Hz).
%      arrw: 1x3 angle rate random walks (rad/s^2/root-Hz).
%       vrw: 1x3 velocity random walks (m/s^2/root-Hz).
%      vrrw: 1x3 velocity rate random walks (m/s^3/root-Hz).
%     g_std: 1x3 gyros standard deviations (radians/s).
%     a_std: 1x3 accrs standard deviations (m/s^2).
%    gb_sta: 1x3 gyros static biases or turn-on biases (radians/s).
%    ab_sta: 1x3 accrs static biases or turn-on biases (m/s^2).
%    gb_dyn: 1x3 gyros dynamic biases or bias instabilities (radians/s).
%    ab_dyn: 1x3 accrs dynamic biases or bias instabilities (m/s^2).
%   gb_corr: 1x3 gyros correlation times (seconds).
%   ab_corr: 1x3 accrs correlation times (seconds).
%    gb_psd: 1x3 gyros dynamic biases PSD (rad/s/root-Hz).
%    ab_psd: 1x3 accrs dynamic biases PSD (m/s^2/root-Hz);
%      freq: 1x1 sampling frequency (Hz).
% ini_align: 1x3 initial attitude at t(1), [roll pitch yaw] (rad).
% ini_align_err: 1x3 initial attitude errors at t(1), [roll pitch yaw] (rad).

ADIS16405.arw      = 2   .* ones(1,3);     % Angle random walks [X Y Z] (deg/root-hour)
ADIS16405.arrw     = zeros(1,3);           % Angle rate random walks [X Y Z] (deg/root-hour/s)
ADIS16405.vrw      = 0.2 .* ones(1,3);     % Velocity random walks [X Y Z] (m/s/root-hour)
ADIS16405.vrrw     = zeros(1,3);           % Velocity rate random walks [X Y Z] (deg/root-hour/s)
ADIS16405.gb_sta   = 3   .* ones(1,3);     % Gyro static biases [X Y Z] (deg/s)
ADIS16405.ab_sta   = 50  .* ones(1,3);     % Acc static biases [X Y Z] (mg)
ADIS16405.gb_dyn   = 0.007 .* ones(1,3);   % Gyro dynamic biases [X Y Z] (deg/s)
ADIS16405.ab_dyn   = 0.2 .* ones(1,3);     % Acc dynamic biases [X Y Z] (mg)
ADIS16405.gb_corr  = 100 .* ones(1,3);     % Gyro correlation times [X Y Z] (seconds)
ADIS16405.ab_corr  = 100 .* ones(1,3);     % Acc correlation times [X Y Z] (seconds)
ADIS16405.freq     = ref.freq;             % IMU operation frequency [X Y Z] (Hz)
% ADIS16405.m_psd     = 0.066 .* ones(1,3);  % Magnetometer noise density [X Y Z] (mgauss/root-Hz)

% ref dataset will be used to simulate IMU sensors.
ADIS16405.t = ref.t;                       % IMU time vector
dt = mean(diff(ADIS16405.t));              % IMU sampling interval

imu1 = imu_si_errors(ADIS16405, dt);       % Transform IMU manufacturer error units to SI units.

imu1.ini_align_err = [0.1 0.1 0.1] .* D2R;                   % Initial attitude align errors for matrix P in Kalman filter, [roll pitch yaw] (radians)  
imu1.ini_align = [ref.roll(1) ref.pitch(1) ref.yaw(1)]; % Initial attitude align at t(1) (radians).

%% ADIS16488 IMU error profile

ADIS16488.arw      = 0.3  .* ones(1,3);     % Angle random walks [X Y Z] (deg/root-hour)
ADIS16488.arrw     = zeros(1,3);            % Angle rate random walks [X Y Z] (deg/root-hour/s)
ADIS16488.vrw      = 0.029.* ones(1,3);     % Velocity random walks [X Y Z] (m/s/root-hour)
ADIS16488.vrrw     = zeros(1,3);            % Velocity rate random walks [X Y Z] (deg/root-hour/s)
ADIS16488.gb_sta   = 0.2  .* ones(1,3);     % Gyro static biases [X Y Z] (deg/s)
ADIS16488.ab_sta   = 16   .* ones(1,3);     % Acc static biases [X Y Z] (mg)
ADIS16488.gb_dyn   = 6.5/3600  .* ones(1,3);% Gyro dynamic biases [X Y Z] (deg/s)
ADIS16488.ab_dyn   = 0.1  .* ones(1,3);     % Acc dynamic biases [X Y Z] (mg)
ADIS16488.gb_corr  = 100  .* ones(1,3);     % Gyro correlation times [X Y Z] (seconds)
ADIS16488.ab_corr  = 100  .* ones(1,3);     % Acc correlation times [X Y Z] (seconds)
ADIS16488.freq     = ref.freq;              % IMU operation frequency [X Y Z] (Hz)
% ADIS16488.m_psd = 0.054 .* ones(1,3);       % Magnetometer noise density [X Y Z] (mgauss/root-Hz)

% ref dataset will be used to simulate IMU sensors.
ADIS16488.t = ref.t;                        % IMU time vector
dt = mean(diff(ADIS16488.t));               % IMU sampling interval

imu2 = imu_si_errors(ADIS16488, dt);        % Transform IMU manufacturer error units to SI units.

imu2.ini_align_err = [1 1 1] .* D2R;                     % Initial attitude align errors for matrix P in Kalman filter, [roll pitch yaw] (radians)  
imu2.ini_align = [ref.roll(1) ref.pitch(1) ref.yaw(1)];  % Initial attitude align at t(1) (radians).

%% dvl: simple model: velocity in NED frame with small additive noise. 

dvl.std= 0.1;
dvl.v= ref.v+0*(dvl.std^2)*randn;
%% IMU1 SYNTHETIC DATA

rng('shuffle')                  % Reset pseudo-random seed

if strcmp(IMU1_DATA, 'ON')      % If simulation of IMU1 data is required ...
    
    fprintf('NaveGo: generating IMU1 ACCR synthetic data... \n')
    
    fb = acc_gen (ref, imu1);   % Generate acc in the body frame
    imu1.fb = ref.fb_syn;
    
    figure
    subplot(3,1,1)
    plot(t_store,imu1.fb(:,1))
    subplot(3,1,2)
    plot(t_store,imu1.fb(:,2))
    subplot(3,1,3)
    plot(t_store,imu1.fb(:,3))
    
    
    fprintf('NaveGo: generating IMU1 GYRO synthetic data... \n')
    
    wb = gyro_gen (ref, imu1);  % Generate gyro in the body frame
    imu1.wb = ref.wb_syn;
    figure
    subplot(3,1,1)
    plot(t_store,imu1.wb(:,1))
    subplot(3,1,2)
    plot(t_store,imu1.wb(:,2))
    subplot(3,1,3)
    plot(t_store,imu1.wb(:,3))
    save imu1.mat imu1
    
    clear wb fb;
    
else
    fprintf('NaveGo: loading IMU1 data... \n')
    
    load imu1.mat
end

%% IMU2 SYNTHETIC DATA

rng('shuffle')					% Reset pseudo-random seed

if strcmp(IMU2_DATA, 'ON')      % If simulation of IMU2 data is required ...
    
    fprintf('NaveGo: generating IMU2 ACCR synthetic data... \n')
    
    fb = acc_gen (ref, imu2);   % Generate acc in the body frame
    imu2.fb = fb;
    
    fprintf('NaveGo: generating IMU2 GYRO synthetic data... \n')
    
    wb = gyro_gen (ref, imu2);  % Generate gyro in the body frame
    imu2.wb = wb;
    
    save imu2.mat imu2
    
    clear wb fb;
    
else
    fprintf('NaveGo: loading IMU2 data... \n')
    
    load imu2.mat
end

%% Print navigation time

to = (ref.t(end) - ref.t(1));

fprintf('Navigation time is%.2f seconds. \n', to)

%% INS/DVL integration using IMU1

if strcmp(IMU1_INS, 'ON')
    
    fprintf('NaveGo: INS/DVL navigation estimates for IMU1... \n')
    
    % Execute INS/DVL integration
    % ---------------------------------------------------------------------
    nav1_e = ins_dvl(imu1, dvl, 'dcm',ref);
    % ---------------------------------------------------------------------
    
    save nav1_e.mat nav1_e
    
else
    
    fprintf('NaveGo: loading INS/GNSS integration for IMU1... \n')
    
    load nav1_e.mat
end


%% Interpolate INS/DVL dataset 

% INS/DVL interpolated according to the reference dataset.

[nav1_r, ref_1] = navego_interpolation (nav1_e, ref);


%% PLOT

figure
plot3(ref.lat,ref.lon,ref.h)
grid minor
zlim([-7,-4])
figure
subplot(3,1,1)
plot(nav1_e.t,ref.roll,nav1_e.t,nav1_e.roll)
xlabel('Time [s]')
ylabel('\phi')

subplot(3,1,2)
plot(nav1_e.t,ref.pitch,nav1_e.t,nav1_e.pitch)
xlabel('Time [s]')
ylabel('\theta')

subplot(3,1,3)
plot(nav1_e.t,ref.yaw,nav1_e.t,nav1_e.yaw)
xlabel('Time [s]')
ylabel('\psi')

figure
subplot(3,1,1)
plot(nav1_e.t,ref.v(:,1),nav1_e.t,nav1_e.vel(:,1))
xlabel('Time [s]')
ylabel('v^N')

subplot(3,1,2)
plot(nav1_e.t,ref.v(:,2),nav1_e.t,nav1_e.vel(:,2))
xlabel('Time [s]')
ylabel('v^E')

subplot(3,1,3)
plot(nav1_e.t,ref.v(:,3),nav1_e.t,nav1_e.vel(:,3))
xlabel('Time [s]')
ylabel('v^D')

pp = abs(nav1_e.Pp(:, 1:16:end).^(0.5)); % Only take diagonal elements from Pp
pi = abs(nav1_e.Pi(:, 1:16:end).^(0.5));
          
     
    % TRAJECTORY
    figure;
    plot3(ref.lon.*R2D, ref.lat.*R2D, ref.h, '--k')
    hold on
    plot3(nav1_e.lon.*R2D, nav1_e.lat.*R2D, nav1_e.h, 'b')
    plot3(ref.lon(1).*R2D, ref.lat(1).*R2D, ref.h(1), 'or', 'MarkerSize', 10, 'LineWidth', 2)
    axis tight
    title('TRAJECTORIES')
    xlabel('Longitude [deg]')
    ylabel('Latitude [deg]')
    zlabel('Altitude [m]')
    view(-25,35)
    legend('TRUE', 'IMU1')
    grid
