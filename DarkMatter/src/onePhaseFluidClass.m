classdef onePhaseFluidClass < fluidClass

    methods

        function obj = onePhaseFluidClass(input, inputTag)
            obj.m       = zeros(input.sim.numpt, 1);
            obj.mdot    = zeros(input.sim.numpt, 1);
            obj.n       = zeros(input.sim.numpt, 1);
            obj.P       = zeros(input.sim.numpt, 1);
            obj.T       = zeros(input.sim.numpt, 1);
            obj.l       = zeros(input.sim.numpt, 1);
            obj.cg      = zeros(input.sim.numpt, 1);
            obj.rho     = zeros(input.sim.numpt, 1);
            obj.tag     = inputTag;

            obj.MW      = input.(inputTag).MW;
            obj.rho     = input.(inputTag).rho;
            obj.name    = input.(inputTag).name;

            obj.r       = input.design.diameter/2;
        end

        function setInitialConditions(obj, input)
            obj.IC.m    = input.m;
            obj.IC.T    = input.T;
            obj.IC.P    = input.P;
            obj.IC.n    = obj.IC.m/obj.MW;

            obj.initializeMass(obj.IC.m);
            obj.initializeTemperature(obj.IC.T);
            obj.initializePressure(obj.IC.P);
            obj.initializeAmount(obj.IC.n);
        end

        function initializeMass(obj, initValue)
            obj.m(1, 1) = initValue;
        end

        function initializeTemperature(obj, initValue)
            obj.T(1, 1) = initValue;
        end

        function initializePressure(obj, initValue)
            obj.P(1, 1) = initValue;
        end

        function initializeAmount(obj, initValue)
            obj.n(1, 1) = initValue;
        end

        function initializeMdot(obj, initValue)
            obj.mdot(1, 1) = initValue;
        end

        function liquidMass = getM(obj, i)
            liquidMass = obj.m(i, 1);
        end

        function liquidPressure = getP(obj, i)
            liquidPressure = obj.P(i, 1);
        end

        function liquidTemperature = getT(obj, i)
            liquidTemperature = obj.T(i, 1);
        end

        function liquidAmount = getN(obj, i)
            liquidAmount = obj.n(i, 1);
        end
    end
end
