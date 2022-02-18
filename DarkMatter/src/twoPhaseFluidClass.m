classdef twoPhaseFluidClass < fluidClass

    methods (Access = public)

        function obj = twoPhaseFluidClass(input, inputTag)
            obj.m       = zeros(input.sim.numpt, 2);
            obj.mdot    = zeros(input.sim.numpt, 2);
            obj.n       = zeros(input.sim.numpt, 2);
            obj.P       = zeros(input.sim.numpt, 2);
            obj.T       = zeros(input.sim.numpt, 2);
            obj.l       = zeros(input.sim.numpt, 2);
            obj.cg      = zeros(input.sim.numpt, 2);
            obj.rho     = zeros(input.sim.numpt, 2);

            obj.MW      = input.(inputTag).MW;
            obj.rho     = input.(inputTag).rho;
            obj.name    = input.(inputTag).name;

            obj.r       = input.design.diameter/2;

            obj.tag     = inputTag;
        end

        function setInitialConditions(obj, input)
            obj.IC.m        = input.m;
            obj.IC.T        = input.T;
            obj.IC.P        = input.P;
            obj.IC.n        = input.m/obj.MW;

            obj.initializeMass(obj.IC.m);
            obj.initializeTemperature(obj.IC.T);
            obj.initializePressure(obj.IC.P);
            obj.initializeAmount(obj.IC.n);
        end

        function initializeMass(obj, initValue)
            obj.m(1, :)     = initValue;
        end

        function initializeTemperature(obj, initValue)
            obj.T(1, :)     = initValue;
        end

        function initializePressure(obj, initValue)
            obj.P(1, :)     = initValue;
        end

        function initializeAmount(obj, initValue)
            obj.n(1, :)     = initValue;
        end

        function initializeMdot(obj, initValue)
            obj.mdot(1, :)  = initValue;
        end

        function liquidMass         = getM(obj, i)
            liquidMass = obj.m(i, :);
        end

        function liquidPressure     = getP(obj, i)
            liquidPressure = obj.P(i, :);
        end

        function liquidTemperature  = getT(obj, i)
            liquidTemperature = obj.T(i, :);
        end

        function liquidAmount       = getN(obj, i)
            liquidAmount = obj.n(i, :);
        end
    end
end
