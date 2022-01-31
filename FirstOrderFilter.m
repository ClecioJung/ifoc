%---------------------------------------------
% Author: Clecio Jung
% First-Order Filter
%---------------------------------------------
classdef FirstOrderFilter < handle
    properties
        % Auxiliary variables
        previousInput
        previousOutput
        % Constants
        zero
        pole
        Ta
        coef
    end
    methods
        % Constructor
        function this = FirstOrderFilter(Ta, zero, pole)
            % Initializes auxiliary variables
            this.previousInput = 0;
            this.previousOutput = 0;
            % Save analog pole and zero values
            this.pole = pole;
            this.zero = zero;
            % Saves sampling time
            this.Ta = Ta;
            % Filter coefficients
            this.coef(1) = (2-pole*this.Ta)/(2+pole*this.Ta);
            this.coef(2) = (pole/zero)*((2+zero*this.Ta)/(2+pole*this.Ta));
            this.coef(3) = (pole/zero)*((zero*this.Ta-2)/(2+pole*this.Ta));
        end
        
        % Calculate control output
        function output = process(this, input)
            this.previousOutput = this.coef(1)*this.previousOutput + ...
                this.coef(2)*input + ...
                this.coef(3)*this.previousInput;
            this.previousInput = input;
            output = this.previousOutput;
        end
    end
end
%---------------------------------------------