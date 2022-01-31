%---------------------------------------------
% Author: Clecio Jung
% PI-D controller with pre-filter
%---------------------------------------------
classdef PI_D < handle
    properties
        % Auxiliary variables
        feedBack
        error
        controlSignal
        % Constants
        xi
        wn
        scl
        minZ
        maxZ
        zero
        Ta
        Kp
        Ti
        Td
        N
        PID
    end
    methods
        % Constructor
        function this = PI_D(Ta)
            % Initializes auxiliary variables
            this.feedBack = zeros(3,1);
            this.error = zeros(3,1);
            this.controlSignal = zeros(3,1);
            % Saves sampling time
            this.Ta = Ta;
        end
        
        % Design controller for desired closed loop characteristics
        function closedLoopResponseProject(this, Mov, ts2, Kma, wma)
            % Save desired closed loop characteristics
            this.xi = abs(log(Mov/100)/sqrt(log(Mov/100)^2 + pi^2));
            this.wn = 4/(this.xi*ts2);
            this.scl = -this.xi*this.wn + 1i*this.wn*sqrt(1 - this.xi^2);
            % Determines the closed-loop zero location
            this.minZ = this.wn/(2*this.xi);
            this.maxZ = this.wn^2/(2*this.xi*this.wn-wma);
            this.zero = sqrt(this.minZ*this.maxZ);
            %this.zero = this.minZ + 0.2*(this.maxZ-this.minZ);
            % Derivative filter coefficient
            this.N = 10;
            % Calculates analog gains
            this.Ti = 1/this.zero;
            this.Kp = this.wn/(Kma*(2*this.xi*this.zero-this.wn));
            this.Td = this.zero/(this.wn^2) - 1/(this.Kp*Kma*wma);
            % Discrete gains
            this.PID(1) = (this.Ta+2*this.Td/this.N)/(this.Ta+this.Td/this.N);
            this.PID(2) = -(this.Td/this.N)/(this.Ta+this.Td/this.N);
            this.PID(3) = this.Kp*(2*this.Ti+this.Ta)/(2*this.Ti);
            this.PID(4) = -this.Kp*(2*this.Ti*this.Ta-this.Ta^2+4*this.Ti*this.Td/this.N)/(2*this.Ti*(this.Ta+this.Td/this.N));
            this.PID(5) = this.Kp*((2*this.Ti-this.Ta)/(2*this.Ti))*((this.Td/this.N)/(this.Ta+this.Td/this.N));
            this.PID(6) = -this.Kp*(this.Td/(this.Ta+this.Td/this.N));
            this.PID(7) = -2*this.PID(6);
            this.PID(8) = this.PID(6);
        end
        
        % Bode plot
        function bode(this)
            PIs_r = tf(this.Kp*[this.Ti 1], [this.Ti 0]);
            PIz_r = tf([this.PID(3) this.PID(4) this.PID(5)], [1 -this.PID(1) -this.PID(2)], this.Ta);
            PIs_y = tf(-this.Kp*[((this.N+1)/this.N)*(this.Ti*this.Td) (this.Ti+this.Td/this.N) 1], this.Ti*[this.Td/this.N 1 0]);
            PIz_y = tf([(this.PID(6)-this.PID(3)) (this.PID(7)-this.PID(4)) (this.PID(8)-this.PID(5))], [1 -this.PID(1) -this.PID(2)], this.Ta);
            figure('WindowState','maximized');
            bode(PIs_r);
            grid on;
            hold on;
            bode(PIz_r);
            bode(PIs_y);
            bode(PIz_y);
            legend('PI(s) (reference)', 'PI(z) (reference)', 'PI(s) (feedback)', 'PI(z) (feedback)');
            set(findall(gcf,'type','line'),'linewidth',2);
        end
        
        % Calculate control output
        function Output = control(this, setPoint, feedBack)
            % Error
            this.feedBack(3) = feedBack;
            this.error(3) = setPoint - this.feedBack(3);
            % Control signal
            this.controlSignal(3) = this.PID(1)*this.controlSignal(2) + ...
                this.PID(2)*this.controlSignal(1) + ...
                this.PID(3)*this.error(3) + ...
                this.PID(4)*this.error(2) + ...
                this.PID(5)*this.error(1) + ...
                this.PID(6)*this.feedBack(3) + ...
                this.PID(7)*this.feedBack(2) + ...
                this.PID(8)*this.feedBack(1);
            % Delay
            this.feedBack(1) = this.feedBack(2);
            this.feedBack(2) =  this.feedBack(3);
            this.error(1) = this.error(2);
            this.error(2) =  this.error(3);
            this.controlSignal(1) = this.controlSignal(2);
            this.controlSignal(2) =  this.controlSignal(3);
            % Output
            Output = this.controlSignal(3);
        end
    end
end
%---------------------------------------------