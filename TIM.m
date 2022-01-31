%---------------------------------------------
% Author: Clecio Jung
% Three-phase induction motor (TIM) in alpha-beta coordinate system
%---------------------------------------------
classdef TIM < handle
    properties
        rko
        n
        m
        fo
        nf
        mf
        % Configuration
        synchronous_referential
        wref
        epslon
        % Runge-Kutta Constants
        ct
        ck
        cx
        % Matrices
        Ad
        Bd
        Af
        Bf
        % Constants
        dt
        p
        gamma
        delta
        Lsigmas
        eta
        Km
        D
        J
        T0
        ba
        % States
        X
        Xf
    end
    
    methods
        function this = TIM(dt,p,gamma,delta,Lsigmas,eta,Km,D,J,T0,ba)
            this.rko = 4;   % Runge-Kutta order
            this.n  = 6;	% Number of states
            this.m  = 5;	% Number of inputs
            
            % Configuration
            this.synchronous_referential = 1;
            this.wref = 0;
            this.epslon = 0.001;
            
            % Saves constants
            this.dt = dt;
            this.p = p;
            this.gamma = gamma;
            this.delta = delta;
            this.Lsigmas = Lsigmas;
            this.eta = eta;
            this.Km = Km;
            this.D = D;
            this.J = J;
            this.T0 = T0;
            this.ba = ba;
            
            switch (this.rko)
                case 1      % Euler Method
                    this.ct = this.dt*0;
                    this.ck = this.dt*0;
                    this.cx = this.dt;
                case 2      % Método de Heun Method
                    this.ct = this.dt*[0; 1];
                    this.ck = this.dt*[0 1; 0 0];
                    this.cx = this.dt*[1/2; 1/2];
                case 3      % Runge-Kutta 3 Method
                    this.ct = this.dt*[0; 1/2; 1];
                    this.ck = this.dt*[0 1/2 -1; 0 0 2; 0 0 0];
                    this.cx = this.dt*[1/6; 4/6; 1/6];
                otherwise   % Runge-Kutta 4 Method
                    this.rko = 4;
                    this.ct = this.dt*[0; 1/2; 1/2; 1];
                    this.ck = this.dt*[0 1/2 0 0; 0 0 1/2 0; 0 0 0 1; 0 0 0 0];
                    this.cx = this.dt*[1/6; 2/6; 2/6; 1/6];
            end
            
            % Matrix Ad
            this.Ad = zeros(this.n,this.n);
            this.Ad(1,1) = -gamma;
            this.Ad(2,2) = this.Ad(1,1);
            this.Ad(3,1) = eta;
            this.Ad(4,2) = this.Ad(3,1);
            this.Ad(1,3) = eta*delta;
            this.Ad(2,4) = this.Ad(1,3);
            this.Ad(3,3) = -eta;
            this.Ad(4,4) = this.Ad(3,3);
            this.Ad(5,5) = -(D/J);
            % Matrix Bd
            this.Bd = zeros(this.n,this.m);
            this.Bd(1,1) = 1/Lsigmas;
            this.Bd(2,2) = this.Bd(1,1);
            this.Bd(5,3) = p/(2*J);
            this.Bd(5,4) = -this.Bd(5,3);
            this.Bd(6,5) = 1;
            
            % Initialize states
            this.X = zeros(this.n,1);
        end
        
        function dX = dX(this,X,vbeta,valpha,Tl)
            % Determines the value of the torque
            Te = this.Km*(X(4)*X(1)-X(3)*X(2));
            Tload = Tl + sign(X(5))*this.T0 + this.ba*X(5)^2;
            
            % Determines the value of the voltages in the rotor flux frame
            cost = cos(X(6));
            sint = sin(X(6));
            vqs = cost*vbeta - sint*valpha;
            vds = sint*vbeta + cost*valpha;
            
            % Referential speed
            if this.synchronous_referential
                if (abs(X(4)) > this.epslon)
                    w = X(5) + this.eta*X(1)/X(4);
                else
                    w = X(5);
                end
            else
                w = this.wref;
            end
            
            % Update matrix Ad
            this.Ad(1,1) = -this.gamma;
            this.Ad(2,2) = this.Ad(1,1);
            this.Ad(2,1) = w;
            this.Ad(1,2) = -this.Ad(2,1);
            this.Ad(2,3) = X(5)*this.delta;
            this.Ad(1,4) = -this.Ad(2,3);
            this.Ad(4,3) = (w-X(5));
            this.Ad(3,4) = -this.Ad(4,3);
            this.Ad(3,1) = this.eta;
            this.Ad(4,2) = this.Ad(3,1);
            this.Ad(1,3) = this.eta*this.delta;
            this.Ad(2,4) = this.Ad(1,3);
            this.Ad(3,3) = -this.eta;
            this.Ad(4,4) = this.Ad(3,3);
            
            % Determines the value of the inputs
            U = [vqs; vds; Te; Tload; w];
            
            % Return values
            dX = this.Ad*X + this.Bd*U;
        end
        
        function [ibeta,ialpha,wr,Te] = process(this,vbeta,valpha,Tl)
            
            % Initializes auxiliary variables for Runge-Kutta method
            kRK = zeros(this.n,this.rko);
            
            for i = 1:this.rko
                % Update the state for Runge-Kutta method
                xRK = this.X + kRK*this.ck(:,i);
                % Determines Runge-Kutta 'K' constants
                kRK(:,i) = this.dX(xRK,vbeta,valpha,Tl);
            end
            
            % Output state
            this.X = this.X + kRK*this.cx;
            
            % Check limits of theta
            if (this.X(6) > 2*pi)
                this.X(6) = this.X(6) - 2*pi;
            elseif (this.X(6) < 0)
                this.X(6) = this.X(6) + 2*pi;
            end
            
            % Determines the value of the currents in the stationary frame
            cost = cos(this.X(6));
            sint = sin(this.X(6));
            ibeta = cost*this.X(1) + sint*this.X(2);
            ialpha = -sint*this.X(1) + cost*this.X(2);
            
            % Return values
            wr = this.X(5);
            Te = this.Km*(this.X(4)*this.X(1)-this.X(3)*this.X(2));
        end
    end
end

