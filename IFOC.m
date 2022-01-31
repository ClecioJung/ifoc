%---------------------------------------------
% Author: Clecio Jung
% IFOC controller
%---------------------------------------------
classdef IFOC < handle
    properties
        % Auxiliary variables
        theta
		previousW		 
        % Objects
        piIqs
        piIds
        piIdm
        Idm
        % Configuration
        use_piIdm
        use_feedforward
        tustin_method
        epslon
        mf
        % Constants
        Ta
        p
        gamma
        delta
        Lsigmas
        eta
        Km
    end
    methods
        % Constructor
        function this = IFOC(Ta,p,gamma,delta,Lsigmas,eta)
            % Initializes auxiliary variables
            this.theta = 0;
			this.previousW = 0;				   
            % Saves constants
            this.Ta = Ta;
            this.p = p;
            this.gamma = gamma;
            this.delta = delta;
            this.Lsigmas = Lsigmas;
            this.eta = eta;
            this.Km = (3*this.p/4)*(delta*Lsigmas);
            
            % Configuration
            this.use_piIdm = 1;
            this.use_feedforward = 0;
			this.tustin_method = 1;					   
            this.epslon = 0.01;
            this.mf = 2; % Relation between wmf and wma
            
            % Current Controllers (IFOC)
            this.piIqs = PI(Ta);
            this.piIqs.poleZeroCancelationProject(this.mf*gamma, 1/(gamma*Lsigmas), gamma);
            this.piIds = PI(Ta);
            this.piIds.poleZeroCancelationProject(this.mf*gamma, 1/(gamma*Lsigmas), gamma);
            % Magnetizing Current Controller
            this.piIdm = PI(Ta);
            this.piIdm.poleZeroCancelationProject(this.mf*eta, 1, eta);
            % Magnetizing Current observer
            this.Idm = Observer(Ta, 1, eta);
        end
        
        function [fq,fd] = park(this,fbeta,falpha)
            ct = cos(this.theta);
            st = sin(this.theta);
            fq = ct*fbeta - st*falpha;
            fd = st*fbeta + ct*falpha;
        end
        
        function [fbeta,falpha] = parkInv(this,fq,fd)
            ct = cos(this.theta);
            st = sin(this.theta);
            fbeta = ct*fq + st*fd;
            falpha = -st*fq + ct*fd;
        end
        
        % Calculate control output
        function [vbeta,valpha,iqs,ids,idm,spIqs,spIds,w] = control(this,spIdm,spTe,wr,ibeta,ialpha)
            % Determines the value of the currents in the rotor flux frame
            [iqs,ids] = this.park(ibeta,ialpha);
            
            % Magnetizing Current observer
            idm = this.Idm.observe(ids);
            
            % Magnetizing Current controller
            if this.use_piIdm
                spIds = this.piIdm.control(spIdm,idm);
            else
                spIds = spIdm;
            end
            
            % Estimates angular position of rotor magnetic flux and determines
            % set point for q axis current
            if (abs(idm) < this.epslon)
                w = wr;
                spIqs = 0;
            else
                w = wr + this.eta*(iqs/idm);
                spIqs = spTe/(this.Km*idm);
            end
			% Integrator - Tustin method
            if this.tustin_method
                this.theta = this.theta + (this.Ta/2)*(w + this.previousW);
                this.previousW = w;
            else % Integrator - Euler method
                this.theta = this.theta + w*this.Ta;
            end
            % Check limits of theta
            if (this.theta > 2*pi)
                this.theta = this.theta - 2*pi;
            elseif (this.theta < 0)
                this.theta = this.theta + 2*pi;
            end
            
            % Current controllers
            uq = this.piIqs.control(spIqs,iqs);
            ud = this.piIds.control(spIds,ids);
            
            % Feedforward compensantion
            if this.use_feedforward
                eq = this.Lsigmas*(w*ids + this.delta*wr*idm);
                ed = -this.Lsigmas*(w*iqs + this.delta*this.eta*idm);
                vqc = uq + eq;
                vdc = ud + ed;
            else
                vqc = uq;
                vdc = ud;
            end
            
            % Return alpha-beta voltages
            [vbeta,valpha] = this.parkInv(vqc,vdc);
        end
    end
end
%---------------------------------------------