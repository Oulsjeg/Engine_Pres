classdef    rocketClass < handle
    properties (Access = public)
        input
        utilities
        propulsion
        airframe
        prop
    end

    methods
        % rocketClass constructor
        function obj = rocketClass(input)
            if (nargin == 0)
                error('ERROR: rocketClass constructor executed without input arguments');
            else
                disp('Message: Building rocketClass');
                obj.input       = input;
                obj.propulsion  = propulsionClass(input);
                obj.airframe    = airframeClass(input);
                obj.utilities   = utilitiesClass(input);

            end
        end
        
        function setMassFlowRates(obj)
            
            for i = 1:obj.propulsion.design.numTanks
                
                switch obj.propulsion.propellants{i, 1}.type
                    case 'Oxidizer'
                        obj.propulsion.performance.Mdotox = obj.propulsion.propellants{i, 1}.propellantTank.propellant.mdot(:, 1);
                    case 'Fuel'
                        obj.propulsion.performance.Mdotf = obj.propulsion.propellants{i, 1}.propellantTank.propellant.mdot(:, 1);
                    otherwise
                        error('rocketClass > setMassFlowRates(): Unknown propellant type.');
                end
                
                obj.propulsion.performance.Mdot = obj.propulsion.performance.Mdotox + obj.propulsion.performance.Mdotf;
            end
        end
    end
end
