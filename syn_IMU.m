function [fb,wb] = syn_IMU(ref.mat)

% INPUT:
%		imu, data structure with IMU error profile in manufacturer units.
%         imu.arw:      angle random walks [X Y Z] (deg/root-hour)
%         imu.arrw:     angle rate random walks [X Y Z] (deg/root-hour/s)
%         imu.vrw:      velocity random walks [X Y Z] (m/s/root-hour)
%		      imu.vrrw:     velocity rate random walks [X Y Z] (deg/root-hour/s)
%         imu.gb_sta:   gyro static biases [X Y Z] (deg/s)
%         imu.ab_sta:   acc static biases [X Y Z] (mg)
%         imu.gb_dyn:   gyro dynamic biases [X Y Z] (deg/s)
%         imu.ab_dynt:  acc dynamic biases [X Y Z] (mg)
%         imu.gb_corr:  gyro correlation times [X Y Z] (seconds)
%         imu.ab_corr:  acc correlation times [X Y Z] (seconds)
%         imu.m_psd:    magnetometer noise density [X Y Z] (mgauss/root-Hz)
%		dt:  IMU sampling interval.
%
% OUTPUT:
%		imu_si: data structure with IMU error profile in SI units.






end

