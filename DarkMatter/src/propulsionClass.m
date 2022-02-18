classdef    propulsionClass < handle

    properties (Access = public)
        info
        input
        design
        utilities
        settings
    end

    properties (Access = public)
        propellants
        performance
    end

    methods
        % propulsionClass constructor
        function obj    = propulsionClass(input)
            if (nargin == 0)            % if no input arguments are passed, initialize as empty arrays
                error('ERROR: propulsionClass constructor executed without input arguments');
            else
                disp('Message: Building propulsionClass');

                obj.input                           = input;
                obj.utilities                       = utilitiesClass(input);

                obj.info.name                       = input.engine.name;
                obj.info.manufacturer               = input.engine.Manu;
                obj.info.engfile                    = input.engine.engfile;
                obj.info.oxidizerTag                = [];
                obj.info.fuelTag                    = [];

                obj.design                          = input.design;
                obj.settings                        = input.settings;

                obj.design.cc.P                     = input.design.Pcc;
                obj.design.cc.R                     = 0;
                obj.design.cc.R                     = 0;
                obj.design.cc.lStar                 = 0;
                obj.design.cc.V                     = 0;
                obj.design.cc.A                     = 0;

                obj.design.nozzle.Rt                = 0;
                obj.design.nozzle.Re                = 0;
                obj.design.nozzle.At                = 0;
                obj.design.nozzle.Ae                = 0;
                obj.design.nozzle.eps               = obj.design.SupArInit;
                obj.design.nozzle.L                 = 0;

                obj.performance.Mdotox              = [];
                obj.performance.Mdotf               = [];
                obj.performance.Mdot                = [];
                obj.performance.Mdotox              = [];
                obj.performance.Pcc                 = [];
                obj.performance.Tcc                 = [];
                obj.performance.cg                  = [];

                obj.getNumTanks(input);
                obj.createPropellantTanks(input);
            end
        end

        function getNumTanks(obj, input)

            inputField      = fieldnames(input);
            obj.design.numTanks        = 0;

            for i = 1:length(inputField)

                fields      = input.(inputField{i, 1});

                if (any(strcmp(fieldnames(fields),'isPropellant')))
                    obj.design.numTanks = obj.design.numTanks + 1;
                end
            end
        end

        function createPropellantTanks(obj, input)

            obj.propellants = cell(obj.design.numTanks, 1);
            tankInfo        = cell(obj.design.numTanks, 1);
            inputField      = fieldnames(input);

            for i = 1:length(inputField)

                inputTag    = inputField{i, 1};
                fields      = input.(inputTag);

                if any(strcmp(fieldnames(fields),'isPropellant')) && fields.('isPropellant')

                    tankInfo{fields.order, 1}   = fields;

                    if (strcmp(tankInfo{fields.order, 1}.fluidtype, 'Oxidizer'))
                        
                        obj.propellants{fields.order, 1} = tankSystemClass(inputTag, input);
                        obj.propellants{fields.order, 1}.order = fields.order;
                        obj.propellants{fields.order, 1}.type  = 'Oxidizer';
                        obj.info.oxidizerTag = inputTag;

                    elseif (strcmp(tankInfo{fields.order, 1}.fluidtype, 'Fuel'))
                        
                        obj.propellants{fields.order, 1} = tankSystemClass(inputTag, input);
                        obj.propellants{fields.order, 1}.order = fields.order;
                        obj.propellants{fields.order, 1}.type  = 'Fuel';
                        obj.info.fuelTag = inputTag;
                        
                    else
                        error('Fluid type %s has not been specialized.', tankInfo{i, 1}.fluidtype);
                    end
                end
            end
        end

        function getCG(obj, data)

            for i = 1:obj.design.numTanks

                obj.propellants{i, 1}.setupTankSystem(data);
                obj.propellants{i, 1}.getCG();
            end
        end

        function output = getCombustion(obj, OF, Pcc, SupAr)
            
            fuel            = {obj.input.(obj.info.fuelTag).name,       ...
                               obj.input.(obj.info.fuelTag).frac,       ...,
                               obj.input.(obj.info.fuelTag).Tinit};

            ox              = {obj.input.(obj.info.oxidizerTag).name,   ...
                               obj.input.(obj.info.oxidizerTag).frac,   ...,
                               obj.input.(obj.info.oxidizerTag).Tinit};

            cea_out                 = CEA_Run(OF, Pcc, SupAr, fuel, ox);

            output.gamm_e           = cea_out.GAMMAs.EXIT1;
            output.cstar            = cea_out.CSTAR.THROAT;
            output.rho_e            = cea_out.RHO.EXIT1;
            output.Ma_e             = cea_out.MachNumber.EXIT1;
            output.Isp              = cea_out.Isp.EXIT1 / obj.utilities.g0;
            output.Tcc              = cea_out.T.CHAMBER;
            output.Te               = cea_out.T.EXIT1;
            output.Pe               = cea_out.P.EXIT1 * 100000;
            output.R_e              = output.Pe / (output.rho_e * output.Te);

        end

        function nozzle_params(obj, Rex, Tcc, mdot)

            Pcc                     = obj.input.design.Pcc * obj.input.settings.cnv;
            gamma_des               = obj.input.design.gamma;
            
            A                       = sqrt(Tcc)/Pcc;
            B                       = sqrt(Rex/gamma_des);
            C                       = ((gamma_des+1)/(2*(gamma_des-1)));
            
            obj.design.nozzle.At    = mdot*A*B*((gamma_des+1)/2)^C;
            obj.design.nozzle.Ae    = obj.design.nozzle.At * obj.design.nozzle.eps;

        end
        
        function Exp_ratio(obj, gamma, P_e, P_0)

            A                       = ((gamma+1)/2)^(1/(gamma-1));
            B                       = (P_e/P_0)^(1/gamma);
            C                       = ((gamma + 1)/(gamma - 1));
            D                       = 1-(P_e/P_0)^((gamma - 1)/gamma);

            obj.design.nozzle.eps   = 1/(A*B*sqrt(C*D));
        end

    end
end
