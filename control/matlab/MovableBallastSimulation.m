close all

phi = 10;  %heel angle degrees
phi_d = 0; 
phi_dd = 0;
gamma = 0;  %arm angle degrees
gamma_d = 0;
gamma_dd = 0;


MBMotorStallTorquekgcm = 100; % kg-cm
MBMotorStallTorque = MBMotorStallTorquekgcm*2.205/2.54; % in-lbf
MBMotorStallCurrent = 28; % Amps
MBMotorFreeCurrent = 5; % Amps
MBMotorFreeSpeed = 86; % rpm
MBmaxAngle = 45;

load = 15; % lbf


stage1 = 16/64;
stage2 = 16/64;
stage1orig = 18/52;
stage2orig = 18/52;

dtime = 0.02;
stopTime = 2;

ka = 120 * 1.2; %torque of keel (N * m)
ks = (120* 1.2^2) + (500 * 0.25); %moment of inertia - drag of keel (rho/2 * area)

ku = 0; %voltage effect on arm acceleration
kf = 0; %frictional resistance to arm acceleration

voltageInput = 0;

phi_array = [];
count = 1;


for t = 0:dtime:stopTime  %simulate for 30 seconds at 0.1 timestep

       
       
       phi_dd = - ka * deg2rad(phi) - ks * phi_d + calcMBRightingMoment(deg2rad(phi), deg2rad(gamma), load);
       gamma_dd = ku * voltageInput - kf * gamma_d - calcMBTorque(deg2rad(phi), deg2rad(gamma), stage1, stage2, load);
 
       phi = phi + dtime * phi_d;  %step forward 
       phi_d = phi_d + dtime * phi_dd;
       gamma = gamma + dtime * gamma_d; 
       gamma_d = gamma_d + dtime*gamma_dd;
       phi_array(count) = phi;
       count = count + 1;
       if gamma > MBmaxAngle
           gamma = MBmaxAngle;
       end
       if gamma < -MBmaxAngle
           gamma = -MBmaxAngle;
       end
       
       
       %animate
       Pm = 5*[sin(deg2rad(phi)),cos(deg2rad(phi))];
       Pb = -3*[sin(deg2rad(gamma)),cos(deg2rad(gamma))];
       
       axis(gca,'equal');
       axis([-3,3,-5,7]);
       
       mast = line([0,Pm(1)],[0,Pm(2)]);
       ballast = line([0,Pb(1)],[0,Pb(2)]);
       
       
       %dispay for time
       pause(dtime)
       
       
       delete(mast);
       delete(ballast);
       
end


close all
plot(0:dtime:stopTime,phi_array);



% tau_m = calcMBTorque(-pi/2, 0, stage1, stage2, load);
%
% pctTauMax = abs(tau_m)/MBMotorStallTorque
%
% maxI = pctTauMax*(MBMotorStallCurrent-MBMotorFreeCurrent)+MBMotorFreeCurrent % Amps
% rpmAtLoad = (1-pctTauMax)*MBMotorFreeSpeed % rpm
% traverseSpeed = rpmAtLoad*stage1*stage2 % rpm


