function [nav_e] = ins_dvl(imu, dvl, att_mode,ref)
% ins_dvl: loosely-coupled integrated navigation system.
% ins_dvl integrates IMU and DVL measurements by using an EKF
%
% INPUT
%   imu: IMU data structure.
%         t: Ix1 time vector (seconds).
%        fb: Ix3 accelerations vector in body frame XYZ (m/s^2).
%        wb: Ix3 turn rates vector in body frame XYZ (radians/s).
%       arw: 1x3 angle random walks (rad/s/root-Hz).
%       vrw: 1x3 velocity random walks (m/s^2/root-Hz).
%      gstd: 1x3 gyros standard deviations (radians/s).
%      astd: 1x3 accrs standard deviations (m/s^2).
%    gb_sta: 1x3 gyros static biases or turn-on biases (radians/s).
%    ab_sta: 1x3 accrs static biases or turn-on biases (m/s^2).
%    gb_dyn: 1x3 gyros dynamic biases or bias instabilities (radians/s).
%    ab_dyn: 1x3 accrs dynamic biases or bias instabilities (m/s^2).
%   gb_corr: 1x3 gyros correlation times (seconds).
%   ab_corr: 1x3 accrs correlation times (seconds).
%    gb_psd: 1x3 gyros dynamic biases PSD (rad/s/root-Hz).
%    ab_psd: 1x3 accrs dynamic biases PSD (m/s^2/root-Hz);
%      freq: 1x1 sampling frequency (Hz).
% ini_align: 1x3 initial attitude at t(1).
% ini_align_err: 1x3 initial attitude errors at t(1).

%	att_mode: attitude mode string.
%      'quaternion': attitude updated in quaternion format. Default value.
%             'dcm': attitude updated in Direct Cosine Matrix format.
%
% OUTPUT
%   nav_e: INS/DVL navigation estimates data structure.
%         t: Ix1 INS time vector (seconds).
%        tg: Mx1 GNSS time vector, when Kalman filter was executed (seconds).
%      roll: Ix1 roll (radians).
%     pitch: Ix1 pitch (radians).
%       yaw: Ix1 yaw (radians).
%       vel: Ix3 NED velocities (m/s).
%       lat: Ix1 latitude (radians).
%       lon: Ix1 longitude (radians).
%         h: Ix1 altitude (m).
%        xi: Mx15 Kalman filter a priori states.
%        xp: Mx15 Kalman filter a posteriori states.
%         z: Mx6  INS/GNSS measurements
%         v: Mx6  Kalman filter innovations.
%         b: Mx6 Kalman filter biases compensations, [gb_dyn ab_dyn].
%         A: Mx225 Kalman filter transition-state matrices, one matrix per
%          row ordered by columns.
%        Pp: Mx225 Kalman filter a posteriori covariance matrices, one
%         matrix per row ordered by columns.
%        Pi: Mx225 Kalman filter a priori covariance matrices, one matrix
%         per row ordered by columns.
%

if nargin < 3, att_mode  = 'quaternion'; end

%% ZUPT detection algorithm

zupt = false;

%% PREALLOCATION

% Constant matrices
I = eye(3);
O = zeros(3);
count=1;
% Length of INS time vector
LI = length(imu.t);

% Length of DVL time vector
LG = length(imu.t);

% Attitude
roll_e  = zeros (LI, 1);
pitch_e = zeros (LI, 1);
yaw_e   = zeros (LI, 1);

% Initialize INS
roll_e(1)  = imu.ini_align(1);
pitch_e(1) = imu.ini_align(2);
yaw_e(1)   = imu.ini_align(3);

DCMnb = euler2dcm([roll_e(1); pitch_e(1); yaw_e(1);]);
DCMbn = DCMnb';
qua   = euler2qua([roll_e(1) pitch_e(1) yaw_e(1)]);

% Velocities
vel_e   = zeros (LI, 3);

% Initialize estimates at INS time = 1
% oracle tells about the initial velocity.
vel_e(1,:) = ref.vel(1,:);

% Positions
lat_e    = zeros (LI,1);
lon_e    = zeros (LI,1);
h_e      = zeros (LI,1);

% Initialize estimates at INS time = 1
h_e(1)   = ref.h(1);
lat_e(1) = ref.lat(1);
lon_e(1) = ref.lon(1);

% Biases
gb_dyn = imu.gb_dyn';
ab_dyn = imu.ab_dyn';

% Initialize Kalman filter matrices

% Prior estimates
kf.xi = [ zeros(1,9), imu.gb_dyn, imu.ab_dyn ]';  % Error vector state
kf.Pi = diag([imu.ini_align_err, dvl.std, dvl.std, dvl.std, 0.01, 0.01, 0.01, imu.gb_dyn, imu.ab_dyn].^2);

kf.R  = 0*diag([dvl.std dvl.std dvl.std].^2);
kf.Q  = diag([imu.arw, imu.vrw, imu.gb_psd, imu.ab_psd].^2);

fb_corrected = (imu.fb(1,:)' + ab_dyn );
fn = (DCMbn * fb_corrected);

% Vector to update matrix F
upd = [dvl.v(1,:) ref.lat(1) ref.h(1) fn'];

% Update matrices F and G
[kf.F, kf.G] = F_update(upd, DCMbn, imu);

[RM,RN] = radius(ref.lat(1));

% Update matrix H
kf.H = [ O I O O O];

kf.z = [dvl.std dvl.std dvl.std]';

% Propagate prior estimates to get xp(1) and Pp(1)
kf = kf_update(kf);


% Kalman filter matrices for later performance analysis
xi = zeros(LG, 15);        % Evolution of Kalman filter a priori states, xi
xp = zeros(LG, 15);        % Evolution of Kalman filter a posteriori states, xp
z = zeros(LG, 3);          % INS/DVL measurements
v = zeros(LG, 3);          % Kalman filter innovations
b  = zeros(LG, 6);         % Biases compensantions after Kalman filter correction

A  = zeros(LG, 225);       % Transition-state matrices, A
Pi = zeros(LG, 225);       % Priori covariance matrices, Pi
Pp = zeros(LG, 225);       % Posteriori covariance matrices, Pp
K  = zeros(LG, 45);       % Kalman gain matrices, K
S  = zeros(LG, 9);       % Innovation matrices, S

% Initialize matrices for Kalman filter performance analysis
xp(1,:) = kf.xp';
Pp(1,:) = reshape(kf.Pp, 1, 225);
b(1,:)  = [imu.gb_sta, imu.ab_sta];
dvlcount = 0;
% INS (IMU) time is the master clock
for i = 2:LI
    
    %% INERTIAL NAVIGATION SYSTEM (INS)
    
    % Print a dot on console every 10,000 INS executions
    
    if (mod(i,1000) == 0), fprintf('. ');  end
    % Print a return on console every 200,000 INS executions
    if (mod(i,2000) == 0), fprintf('\n'); end
    
    % IMU sampling interval
    dti = imu.t(i) - imu.t(i-1);
    
    % Inertial sensors corrected with KF biases estimation
    wb_corrected = (imu.wb(i,:)' + gb_dyn );
    fb_corrected = (imu.fb(i,:)' + ab_dyn );
    
    % Turn-rates update
    omega_ie_n = earthrate(lat_e(i-1));
    omega_en_n = transportrate(lat_e(i-1), vel_e(i-1,1), vel_e(i-1,2), h_e(i-1));
    
    % Attitude update
    [qua_n, DCMbn, euler] = att_update(wb_corrected, DCMbn, qua, ...
        omega_ie_n, omega_en_n, dti, att_mode);
    roll_e(i) = euler(1);
    pitch_e(i)= euler(2);
    yaw_e(i)  = euler(3);
    qua = qua_n;
    
    % Gravity update
    gn = gravity(lat_e(i-1), h_e(i-1));
    
    % Velocity update
    fn = (DCMbn * fb_corrected);
    vel_n = vel_update(fn, vel_e(i-1,:), omega_ie_n, omega_en_n, gn', dti);
    vel_e (i,:) = vel_n;
    
    % Position update
    pos = pos_update([lat_e(i-1) lon_e(i-1) h_e(i-1)], vel_e(i,:), dti);
    lat_e(i) = pos(1);
    lon_e(i) = pos(2);
    h_e(i)   = pos(3);
    
    %% KALMAN FILTER UPDATE
    dtDVL=50;
    dvl.noise=[1 1 1];
    
    dvl.vel = ref.vel(i,:)+(randn/sqrt(dtDVL));
    if (mod(i,50)==0)  && i<40000      
        %% INNOVATIONS
        zv = (vel_e(i,:) - dvl.v(i,:))';
        count=count+1;
        %% KALMAN FILTER
        
        % Vector to update matrix F
        upd = [vel_e(i,:) lat_e(i) h_e(i) fn'];
        
        % Update matrices F and G
        [kf.F, kf.G] = F_update(upd, DCMbn, imu);
        
        % Update matrix H
        
        kf.H = [ O I O O O];
        kf.R = diag((dvl.noise).^2);
        kf.z = zv;
        
        % Execute the extended Kalman filter
        kf.xp(1:9) = 0;              % states 1:9 are forced to be zero (error-state approach)
        kf = kalman(kf, dtDVL);
        
        %% INS/GNSS CORRECTIONS
        
        %         % Quaternion corrections
        %         % Crassidis. Eq. 7.34 and A.174a.
        antm = [0 qua_n(3) -qua_n(2); -qua_n(3) 0 qua_n(1); qua_n(2) -qua_n(1) 0];
        qua = qua_n + 0.5 .* [qua_n(4)*eye(3) + antm; -1.*[qua_n(1) qua_n(2) qua_n(3)]] * kf.xp(1:3);
        qua = qua / norm(qua);       % Brute-force normalization
              
        % Attitude corrections
        roll_e(i)  = roll_e(i)  - kf.xp(1);
        pitch_e(i) = pitch_e(i) - kf.xp(2);
        yaw_e(i)   = yaw_e(i)   - kf.xp(3);
        
        % Velocity corrections
        vel_e (i,1) = vel_e (i,1) - kf.xp(4);
        vel_e (i,2) = vel_e (i,2) - kf.xp(5);
        vel_e (i,3) = vel_e (i,3) - kf.xp(6);
        
        % Position corrections
        lat_e(i) = lat_e(i) - (kf.xp(7));
        lon_e(i) = lon_e(i) - (kf.xp(8));
        h_e(i)   = h_e(i)   - kf.xp(9);
        
        % Biases corrections
        gb_dyn   = kf.xp(10:12);
        ab_dyn   = kf.xp(13:15);
        
        % Matrices for later Kalman filter performance analysis
        xi(i,:) = kf.xi';
        xp(i,:) = kf.xp';
        b(i,:) = [gb_dyn', ab_dyn'];
        A(i,:)  = reshape(kf.A, 1, 225);
        Pi(i,:) = reshape(kf.Pi, 1, 225);
        Pp(i,:) = reshape(kf.Pp, 1, 225);
       
        if(zupt == false)
            v(i,:)  = kf.v';
            K(i,:)  = reshape(kf.K, 1, 45);
            S(i,:)  = reshape(kf.S, 1, 9);
        else
            zupt = false;
            v(i,:)  = [ kf.v' 0 0 0 ]';
            K(i,1:45)  = reshape(kf.K, 1, 45);
            S(i,1:9)  = reshape(kf.S, 1, 9);
        end
    end
end

%% Summary from INS/GNSS integration

nav_e.t     = imu.t(1:i, :);    % INS time vector
nav_e.tg    = imu.t(1:i, :);           % GNSS time vector, which is the time vector when the Kalman filter was executed
nav_e.roll  = roll_e(1:i, :);   % Roll
nav_e.pitch = pitch_e(1:i, :);  % Pitch
nav_e.yaw   = yaw_e(1:i, :);    % Yaw
nav_e.vel   = vel_e(1:i, :);    % NED velocities
nav_e.lat   = lat_e(1:i, :);    % Latitude
nav_e.lon   = lon_e(1:i, :);    % Longitude
nav_e.h     = h_e(1:i, :);      % Altitude

nav_e.xi    = xi;       % A priori states
nav_e.xp    = xp;       % A posteriori states
nav_e.m     = z;        % INS/GNSS measurements
nav_e.v     = v;        % Kalman filter innovations
nav_e.b     = b;        % Biases compensations

nav_e.A     = A;        % Transition matrices
nav_e.Pi    = Pi;       % A priori covariance matrices
nav_e.Pp    = Pp;       % A posteriori covariance matrices
nav_e.K     = K;        % Kalman gain matrices
nav_e.S     = S;        % Innovation matrices

fprintf('\n');

end