%---------------------------------------------
% Author: Clecio Jung
% PI controller
%---------------------------------------------
classdef PI < handle
    properties
        % Auxiliary variables
        controlSignal
        previousError
        % Constants
        wmf
        Ta
        xi
        wn
        scl
        zero
        Kp
        Ti
        Kc
        zc
    end
    methods
        % Constructor
        function this = PI(Ta)
            % Initializes auxiliary variables
            this.controlSignal = 0;
            this.previousError = 0;
            % Saves sampling time
            this.Ta = Ta;
        end
        
        % Design controller for desired closed loop characteristics
        function poleZeroCancelationProject(this, wmf, Kma, wma)
            % Save desired closed loop bandwidth
            this.wmf = wmf;
            % Calculates analog gains using pole zero cancellation
            this.Ti = 1/wma;
            this.Kp = (wmf*this.Ti)/Kma;
            % Discrete gains
            this.discreet(this.Kp, this.Ti);
        end
        
        % Design controller for desired closed loop characteristics
        function closedLoopResponseProject(this, Mov, ts2, Kma, wma)
            % Save desired closed loop characteristics
            this.xi = abs(log(Mov/100)/sqrt(log(Mov/100)^2 + pi^2));
            this.wn = 4/(this.xi*ts2);
            this.scl = -this.xi*this.wn + 1i*this.wn*sqrt(1 - this.xi^2);
            % Calculates analog gains
            this.Kp = (2*this.xi*this.wn - wma)/(Kma*wma);
            this.Ti = (2*this.xi*this.wn - wma)/(this.wn^2);
            % Closed loop zero
            this.zero = 1/this.Ti;
            % Discrete gains
            this.discreet(this.Kp, this.Ti);
        end
        
        % Bode plot
        function bode(this)
            PIs = tf(this.Kp*[this.Ti 1], [this.Ti 0]);
            PIz = tf(this.Kc*[1 -this.zc], [1 -1]);
            figure('WindowState','maximized');
            bode(PIs);
            grid on;
            hold on;
            bode(PIz);
            legend('PI(s)', 'PI(z)');
            set(findall(gcf,'type','line'),'linewidth',2);
        end
        
        % Design controller for desired closed loop characteristics
        function KiSulProject(this, wcl, J, Kol)
            % Save desired closed loop bandwidth
            this.wcl = wcl;
            % Calculates analog gains using Ki-Sul suggestions
            % book Control of Electric Machine Drive Systems,
            % author Seung-Ki Sul, year 2011, cap 4 pag. 200
            this.discreet((wcl*J)/Kol, 10/wcl);
        end
        
        % Discrete gains (Tustin approximation)
        function discreet(this,Kp,Ti)
            this.Kp = Kp;
            this.Ti = Ti;
            this.Kc = this.Kp*(2*this.Ti+this.Ta)/(2*this.Ti);
            this.zc = (2*this.Ti-this.Ta)/(2*this.Ti+this.Ta);
        end
        
        % Calculate control output
        function Output = control(this, setPoint, feedBack)
            error = setPoint - feedBack;
            this.controlSignal = this.controlSignal + this.Kc*(error - this.zc*this.previousError);
            this.previousError = error;
            % Output
            Output = this.controlSignal;
        end
    end
end
%---------------------------------------------