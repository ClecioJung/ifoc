%---------------------------------------------
% Author: Clecio Jung
% 1st order observer without zero (used to estimate motor's flux)
%---------------------------------------------
classdef Observer < handle
    properties
        % Auxiliary variables
        previousInput
        previousOutput
        % Constants
        Kma
        wma
        Ta
        Kd
        zd
    end
    methods
        % Constructor
        function this = Observer(Ta, Kma, wma)
            % Initializes auxiliary variables
            this.previousInput = 0;
            this.previousOutput = 0;
            % Saves sampling time
            this.Ta = Ta;
            % Run observer project
            this.project(Kma, wma);
        end
        
        % Observer design (ZOH approach)
        function project(this, Kma, wma)
            this.Kma = Kma;
            this.wma = wma;
            this.zd = exp(-this.Ta*this.wma);
            this.Kd = this.Kma*(1-this.zd);
        end
        
        % 1st order observer
        function output = observe(this, input)
            output = this.zd*this.previousOutput + this.Kd*(this.previousInput);
            this.previousInput = input;
            this.previousOutput = output;
        end
    end
end
%---------------------------------------------