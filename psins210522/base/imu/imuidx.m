function imu = imuidx(data, idx, gunit, aunit, ts)
% Extract IMU(Gyro/Acc/T) from a data array.
%
% Prototype: imu = imuidx(data, idx, gunit, aunit, ts)
% Inputs: data - data contains IMU
%         idx - column index of IMU/T
%         gunit,aunit - gyro- & acc- unit
%         ts - sampling interval
% Output: imu - =[gyro,acc,t] angular- & velocity- increment & time tag
%
% See also: imurfu, gpsidx, avpidx.

% Copyright(c) 2009-2021, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 31/01/2020
global glv
    imu = data(:,abs(idx));
    if nargin<5, ts=1; end
    if nargin<4, aunit=glv.g0; end
    if nargin<3, gunit=glv.dps; end
    if length(idx)<7, imu(:,7)=(1:length(imu))'*ts; end
    imu(:,1:3) = imu(:,1:3)*gunit*ts;
    imu(:,4:6) = imu(:,4:6)*aunit*ts;
    for k=1:6
        if(idx(k))<0, imu(:,k)=-imu(:,k); end
    end
    