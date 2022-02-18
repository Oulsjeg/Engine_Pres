classdef (Abstract) fluidClass < handle

    properties
        m               % Mass                                      (double array)
        n               % Amount                                    (double array)
        T               % temperature                               (double array)
        P               % Pressure                                  (double array)
        l               % Length of fluid                           (double array)
        mdot            % Mass flow rate                            (double array)
        rho             % Density                                   (double array)
        cg              % Center of gravity                         (double array)
        MW              % Molar mass                                (double)
        r               % Radius                                    (double)
        IC              % Initial conditions                        (struct)
        tag             % Input field name                          (str)
        name            % Fluid name                                (str)

        % Blowdown characteristics handle                           (handle)
        blowdownCharacteristics
    end

    methods

        setInitialConditions(obj)
        initializeMass(obj)
        initializeTemperature(obj)
        initializePressure(obj)
        initializeAmount(obj)
    end
end
