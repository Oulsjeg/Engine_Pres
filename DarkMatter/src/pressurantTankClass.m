classdef pressurantTankClass < handle

    properties

        name            % Propellant name                           (str)
        type            % Propellant type                           (str)
        inputTag        % Tag for name of dedicated input field     (str)
        input           % Input structure                           (struct)
        tank            % Tank body                                 (struct)
        pressurant      % Liquid phase                              (fluidClass)
        cg              % CG of propellant tank system              (double)
        m               % Mass of propellant tank system            (double)
        offset          % Distance offset before tank system        (double)
        utilities       % Utilities                                 (utilitiesClass)
    end

    methods
        function obj        = pressurantTankClass(inputTag, input)

            obj.tank.m          = input.(inputTag).mTank;
            obj.tank.l          = input.(inputTag).lTank;
            obj.tank.cg         = obj.tank.l/2;

            obj.input           = input;
            obj.name            = input.(inputTag).name;
            obj.cg              = zeros(input.sim.numpt, 1);
            obj.offset          = input.(inputTag).offset;
            obj.m               = input.(inputTag).mInit*ones(input.sim.numpt, 1);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Create propellant liquid, vapor and pressurant vapor classes

            % Initialize propellant liquid phase
            obj.pressurant      = onePhaseFluidClass(input, inputTag);

            initStruct.m        = input.(inputTag).mInit;
            initStruct.T        = input.(inputTag).Tinit;
            initStruct.P        = input.(inputTag).Pinit;

            obj.pressurant.setInitialConditions(initStruct);

            obj.utilities       = utilitiesClass(input);
        end

        function setupTankSystem(obj, mdot)

            obj.liq.mdot     = mdot;

            for i = 2:length(obj.pres.mdot)
                obj.liq.m(i) = obj.liq.m(i-1) - obj.input.sim.dt*obj.liq.mdot(i);
            end

            obj.m               = obj.tank.m.*ones(length(obj.liq.m), 1) ...
                                + obj.pres.m;
        end

        function getCG(obj)
            obj.cg(:, 1) = 0.5*obj.tank.l;
        end

    end

end
