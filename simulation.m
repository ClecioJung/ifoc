%---------------------------------------------
% Author: Clecio Jung
% Date: 05/07/2020
%---------------------------------------------

clc;
clear all;
close all;

timerVal = tic();

% Import motor parameters
run('parameters.m');

%---------------------------------------------
% Simulation parameters
%---------------------------------------------

% Load parameters
simParams.Tload = 10;           % Load torque (Nm)
simParams.TlOn  = 5.0;          % Instant time to trigger the load
simParams.TlOff = 10.0;         % Time to deactivate the load

% Simulation parameters
simParams.tf = 15;              % Simulation end time (s)

%---------------------------------------------
% Control parameters
%---------------------------------------------

% Controller options
controlParams.pidSpeed = 1;

% Controller parameters
controlParams.fs = 10000;               % Sampling frequency
controlParams.Ta = 1/controlParams.fs;  % Sampling time (s)
controlParams.wOn = 1.0;                % Instant when speed setpoint is activated (step)
controlParams.wSP = 1.0*params.wn;      % Speed setpoint value (step)
controlParams.spTe = 0.5;        

% Estimated motor parameters
controlParams.p = params.p;
controlParams.gamma = params.gamma;
controlParams.delta = params.delta;
controlParams.Lsigmas = params.Lsigmas;
controlParams.eta = params.eta;
controlParams.D = params.D;
controlParams.J = params.J;
controlParams.Imr = params.Imr;

%---------------------------------------------
% Motor simulation controlled by IFOC algorithm
%---------------------------------------------

% Simulation coefficients
N  = round(simParams.tf/controlParams.Ta);  % Amount of samples for controller
Na = 10;                                    % Number of points between two samples
dt = controlParams.Ta/Na;                   % Simulation step (s)

% TIM object
tim = TIM(dt,params.p,params.gamma,params.delta,params.Lsigmas,params.eta,params.Km,params.D,params.J,params.T0,params.ba);

% Speed Controller
if controlParams.pidSpeed
    pidW = PI_D(controlParams.Ta);
    pidW.closedLoopResponseProject(1.0, 2.0, controlParams.p/(2*controlParams.D), controlParams.D/controlParams.J);
    filterW = FirstOrderFilter(controlParams.Ta, 10*abs(real(pidW.scl)), pidW.zero);
else
    piW = PI(controlParams.Ta);
    piW.closedLoopResponseProject(1.0, 4.0, controlParams.p/(2*controlParams.D), controlParams.D/controlParams.J);
end

% IFOC controller
ifoc = IFOC(controlParams.Ta,controlParams.p,controlParams.gamma,controlParams.delta,controlParams.Lsigmas,controlParams.eta);

% Initialization
data.t = (1:N)*controlParams.Ta;
data.ibeta = zeros(1,N);
data.ialpha = zeros(1,N);
data.iqs = zeros(1,N);
data.ids = zeros(1,N);
data.idm = zeros(1,N);
data.wr = zeros(1,N);
data.vbeta = zeros(1,N);
data.valpha = zeros(1,N);
data.vqs = zeros(1,N);
data.vds = zeros(1,N);
data.vqc = zeros(1,N);
data.vdc = zeros(1,N);
data.Te = zeros(1,N);
data.Tl = zeros(1,N);
data.spIqs = zeros(1,N);
data.spIds = zeros(1,N);
data.spIdm = zeros(1,N);
data.spTe = zeros(1,N);
data.spW = zeros(1,N);
data.spWf = zeros(1,N);
data.uqs = zeros(1,N);

% Initialization
wr = 0;
ibeta = 0;
ialpha = 0;
vbeta = 0;
valpha = 0;
spWf = 0;
motorVbeta = 0;
motorValpha = 0;

% IFOC simulation
for n = 1:N
    time = data.t(n);
    
    % Load torque
    Tl = simParams.Tload*(heaviside(time-simParams.TlOn)-heaviside(time-simParams.TlOff));
    
    % Determines the alpha-beta referenthial values
    [vqs,vds] = ifoc.park(vbeta,valpha);
    
    % Speed controller
    spW = controlParams.wSP*heaviside(time - controlParams.wOn);
    if controlParams.pidSpeed
        spWf = filterW.process(spW);
        spTe = pidW.control(spWf,wr);
    else
        spTe = piW.control(spW,wr);
    end
    
    % IFOC control
    spIdm = controlParams.Imr;
    [vbeta,valpha,iqs,ids,idm,spIqs,spIds,w] = ifoc.control(spIdm,spTe,wr,ibeta,ialpha);
    
    % TIM simulation
    for k = 1:Na
        % Update beta-alpha voltage applyed to the motor
        if (k > floor(Na/2))
            motorVbeta = vbeta;
            motorValpha = valpha;
        end
        
        % Motor + anti-aliasing filter model
        [ibeta,ialpha,wr,Te] = tim.process(motorVbeta,motorValpha,Tl);
    end
    
    % Save simulation values
    data.ibeta(n) = ibeta;
    data.ialpha(n) = ialpha;
    data.iqs(n) = iqs;
    data.ids(n) = ids;
    data.idm(n) = idm;
    data.wr(n) = wr;
    data.w(n) = w;
    data.vbeta(n) = vbeta;
    data.valpha(n) = valpha;
    data.vqs(n) = vqs;
    data.vds(n) = vds;
    data.Te(n) = Te;
    data.Tl(n) = Tl;
    data.spIqs(n) = spIqs;
    data.spIds(n) = spIds;
    data.spIdm(n) = spIdm;
    data.spTe(n) = spTe;
    data.spW(n) = spW;
    data.spWf(n) = spWf;
end

% Power and their relationships
data.P = 3/2*(data.vbeta.*data.ibeta + data.valpha.*data.ialpha);
data.Q = 3/2*(data.vbeta.*data.ialpha - data.valpha.*data.ibeta);
data.Pe = (2/params.p)*data.Te.*data.wr;
data.S = sqrt(data.P.^2+data.Q.^2);
data.fp = min(max(data.P./data.S,0.0),1.0);
data.efficiency = min(max(data.Pe./data.P,0.0),1.0);
data.fp(data.t < controlParams.wOn) = 0.0;
data.efficiency(data.t < controlParams.wOn) = 0.0;
data.Energy = cumsum(data.P.*controlParams.Ta);

% RMS values
data.vrms = sqrt((data.vqs.^2+data.vds.^2)./2);
data.irms = sqrt((data.iqs.^2+data.ids.^2)./2);

%---------------------------------------------
% Charts
%---------------------------------------------

figure;
plot(data.t,data.vqs);
hold on;
grid on;
plot(data.t,data.vds);
xlabel('Time (s)','interpreter','latex');
ylabel('Voltage (V)','interpreter','latex');
title('Controller Voltages','interpreter','latex');
legend('$v_{qs}$','$v_{ds}$','interpreter','latex','Location','east');

figure;
plot(data.t,data.vbeta);
hold on;
grid on;
plot(data.t,data.valpha);
xlabel('Time (s)','interpreter','latex');
ylabel('Voltage (V)','interpreter','latex');
title('Controller Voltages','interpreter','latex');
legend('$v_{\beta s}$','$v_{\alpha s}$','interpreter','latex','Location','east');

figure;
plot(data.t,data.iqs);
hold on;
grid on;
plot(data.t,data.ids);
plot(data.t,data.idm,'--');
xlabel('Time (s)','interpreter','latex');
ylabel('Current (A)','interpreter','latex');
title('Controller Currents','interpreter','latex');
legend('$i_{qs}$','$i_{ds}$','$\widehat{i}_{md}$','interpreter','latex','Location','northeast');

figure;
plot(data.t,data.ibeta);
hold on;
grid on;
plot(data.t,data.ialpha);
xlabel('Time (s)','interpreter','latex');
ylabel('Current (A)','interpreter','latex');
title('Controller Currents','interpreter','latex');
legend('$i_{\beta s}$','$i_{\alpha s}$','interpreter','latex','Location','northeast');

figure;
plot(data.t,data.Te);
hold on;
grid on;
plot(data.t,data.spTe);
plot(data.t,data.Tl);
xlabel('Time (s)','interpreter','latex');
ylabel('Torque (Nm)','interpreter','latex');
title('Torque','interpreter','latex');
legend('$T_e$','$T_e^*$','$T_L$','interpreter','latex','Location','northeast');

figure;
plot(data.t,data.wr);
hold on;
grid on;
plot(data.t,data.spW);
plot(data.t,data.w);
xlabel('Time (s)','interpreter','latex');
ylabel('Speed (rad/s)','interpreter','latex');
title('Speed Controller','interpreter','latex');
legend('$\omega_r$','$\omega_r^*$','$\omega$','interpreter','latex','Location','southeast');

figure;
plot(data.t,data.P);
hold on;
grid on;
plot(data.t,data.Q);
plot(data.t,data.S);
xlabel('Time (s)','interpreter','latex');
ylabel('Power (W)','interpreter','latex');
title('Power','interpreter','latex');
legend('P','Q','S','interpreter','latex','Location','east');

elapsedTime = toc(timerVal);
fprintf('Elapsed Time: %.2f s\n',elapsedTime);

%---------------------------------------------