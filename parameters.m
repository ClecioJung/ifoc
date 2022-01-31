%---------------------------------------------
% Motor parameters (obtained form a 4cv 8-pole three-phase induction motor)
%---------------------------------------------

% Nominal values
params.Vn  = 220;                       % Rated Voltage (V)
params.Pn  = 4*735.499;                 % Rated power (W)
params.In  = 12.7;                      % Rated Current (A)
params.fn  = 60;                        % Rated Frequency (Hz)
params.wn  = 2*pi()*params.fn;          % Angular Frequency (rad/s)
params.Tn  = 33.34;                     % Rated Torque (Nm)
params.p   = 8;                         % Pole number
params.N   = 860;                       % Rated Speed (RPM)
params.Ns  = 120*params.fn/params.p;    % Synchronous Speed (RPM)
params.sn  = 100*(params.Ns-params.N)/params.Ns;    % Rated Slip (%)
% Resistances
params.rs  = 0.467;
params.rr  = 0.380;
% Reactances
params.Xls = 1.244;
params.Xlr = params.Xls;
params.Xm  = 14.957;
params.Xss = params.Xm + params.Xls;
params.Xrr = params.Xm + params.Xlr;
% Mechanical parameters
params.J   = 2*0.1033 + 40.5e-6;
params.D   = 0.0015;
params.T0  = 0.3777;
params.ba  = 3e-7;
% Inductances
params.Lls = params.Xls/params.wn;
params.Llr = params.Xlr/params.wn;
params.Lm = params.Xm/params.wn;
params.Lss = params.Xss/params.wn;
params.Lrr = params.Xrr/params.wn;
% Secondary coefficients
params.sigma = 1 - params.Lm^2/(params.Lss*params.Lrr);
params.delta = params.Lm^2/(params.sigma*params.Lss*params.Lrr);
params.Lsigmas = params.sigma*params.Lss;
params.eta = params.rr/params.Lrr;
params.gamma = (params.rs+(params.Lm/params.Lrr)^2*params.rr)/params.Lsigmas;
params.wmec = params.D/params.J;
% Nominal flux
params.fluxn = sqrt(2)*params.Vn/params.wn;
% Magnetization current
params.Imr = params.fluxn/params.Lm;
% Torque-current proportionality constant
params.Km = (3*params.p/4)*(params.Lm^2/params.Lrr);

%---------------------------------------------