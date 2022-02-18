classdef simulationClass < handle

    properties (Access = public)
        input
        utilities
        rocket
        flight
        t, tBurn, numpt, relax, altConvCrit
    end

    methods
        % simulationClass constructor
        function obj	= simulationClass(input)
            if (nargin == 0)
                error('ERROR: simulationClass constructor executed without input arguments');

            else
                disp('Message: Building simulationClass');
                obj.rocket          = rocketClass(input);

                obj.input           = input;
                obj.utilities       = utilitiesClass(input);

                obj.numpt           = input.sim.numpt;
                obj.tBurn           = input.sim.tBurn;
                obj.relax           = input.sim.relax;
                obj.altConvCrit     = input.sim.altConvCrit;
                obj.t               = linspace(0, obj.tBurn, obj.numpt);

                obj.flight.x        = obj.utilities.zeroArray;
                obj.flight.y        = obj.utilities.zeroArray;
                obj.flight.u        = obj.utilities.zeroArray;
                obj.flight.v        = obj.utilities.zeroArray;
                obj.flight.ax       = obj.utilities.zeroArray;
                obj.flight.ay       = obj.utilities.zeroArray;
                obj.flight.thetaL   = input.sim.thetaL;
                obj.flight.Ma       = 0;
                obj.flight.g        = 0;
                obj.flight.type     = input.settings.flightType;

                obj.setFlightDynamicsType();
            end
        end

        function setFlightDynamicsType(obj)
            switch obj.flight.type
                case '2DOF'
                    obj.flight.get = @ obj.get2DOFflightDynamics;
                otherwise
                    error('Flight dynamics mode %s has not been specialized.', obj.flight.type);
            end
        end

        function get2DOFflightDynamics(obj)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % temporary:
            % 1) use a fake OF value array until momo merges code
            % OF = linspace(3, 3.2, numpt);
            % 2) use a fake fuel temp array fuelTemp = 298*ones(numpt, 1);
            % 3) use a fake ox temp array   oxTemp = 278*ones(numpt, 1);
            % 4) use fake mdot              mdot = obj.design.mdotox * (1 + 1/obj.design.OF)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            prop                    = obj.rocket.propulsion;

            num                     = obj.input.sim.numpt;
            num_OF                  = obj.input.settings.num_OF;

            m_wet                   = 70 * ones(num ,1);
            m_dot                   = prop.performance.Mdot;
            mtot                    = m_wet - m_dot;

            fakeOF                  = linspace(3, 3.2, num)';

            obj.utilities.getTcurve();
            T_curve                 = obj.utilities.T_curve;

            R_air                   = obj.utilities.airR;   %Specific gas constant for dry air (J/(kg.k))
            g                       = -obj.utilities.g0;
            Ar                      = pi * (obj.input.design.diameter/2)^2;  %Drag reference area [m^2]
            eff                     = obj.input.settings.efficiency;

            % Initializing Result Arrays
            thrust                  = zeros(num,1);                 % come back and fix this
            Pcc                     = zeros(num,1);
            Tcc                     = zeros(num,1);
            Isp                     = zeros(num, 1);
            Pe                      = zeros(num, 1);
            m_f                     = zeros(num,1);
            drag                    = zeros(num,1);
            rho_air                 = zeros(num,1);
            velx                    = zeros(num,1);
            vely                    = zeros(num,1);
            vel                     = zeros(num,1);
            sy                      = zeros(num,1);
            sx                      = zeros(num,1);
            Ma                      = zeros(num,1);
            Ay                      = zeros(num,1);
            Ax                      = zeros(num,1);
            gforce                  = zeros(num,1);
            force                   = zeros(num,1);
            cstar_temp              = zeros(num_OF, 1);

            % Initialize drag model
            cd_curve                = obj.rocket.airframe.aero.cd;

            OF_Vec                  = linspace(obj.input.settings.OF_i,obj.input.settings.OF_f, obj.input.settings.num_OF)';
            Pcc_design              = obj.input.design.Pcc * obj.input.settings.cnv;
            SupAr                   = obj.input.design.SupArInit;

            combustion              = prop.getCombustion(obj.input.design.OF, Pcc_design, SupAr);
            Rex                     = combustion.Pe / (combustion.rho_e * combustion.Te);
            mdot                    = prop.design.mDotOx * (1 + 1/prop.design.OF);

            prop.nozzle_params(Rex, combustion.Tcc, mdot);

            A_t                     = prop.design.nozzle.At;

            for i = 1:num_OF

                combustion          = prop.getCombustion(OF_Vec(i), Pcc_design, SupAr);
                cstar_temp(i, 1)    = combustion.cstar;
            end

            cstar_curve             = fit(OF_Vec, cstar_temp, 'cubicinterp');

            Pcc(1)                      = cstar_curve(fakeOF(1)) * m_dot(1) / A_t;
            combustion                  = prop.getCombustion(fakeOF(1), Pcc(1), SupAr);

            alt_corr                    = (combustion.Pe - obj.utilities.getAmbientPressure(sy(1)))*prop.performance.A_e;

            thrust(1,1)                 = prop.performance.Mdot(1,1) * eff * combustion.Isp* 9.81 + alt_corr;

            Ay(1,1)                     = thrust(1,1)*sind(obj.input.sim.thetaL)/ mtot(1,1);
            Ax(1,1)                     = thrust(1,1)*cosd(obj.input.sim.thetaL)/ mtot(1,1);

            % gamma   = gamma_curve(of design);
            % inputs.SupSr = somewhere in design;

            for i = 2:num

                Pcc(i)                  = cstar_curve(fakeOF(i)) * m_dot(i) / A_t;

                % inputs.fuelTemp = fuelTemp(i);
                % inputs.oxTemp   = oxTemp(i);

                combustion              = prop.getCombustion(fakeOF(i), Pcc(i), SupAr);

                Pe(i, 1)                = combustion.Pe;
                Isp(i, 1)               = combustion.Isp;

                m_f(i,1)                = m_f(i-1,1) - prop.performance.Mdot(i-1,1)*obj.input.sim.dt;

                rho_air(i,1)            = obj.utilities.getAmbientPressure(sy(i-1,1)) / (R_air * T_curve(sy(i-1,1)));       % Air density

                alt_corr                = (Pe(i-1) - obj.utilities.getAmbientPressure(sy(i-1,1)) )*prop.performance.A_e;    % Altitude correction

                thrust(i,1)             = prop.performance.Mdot(i-1) * eff*Isp(i-1) * 9.81 + alt_corr;                      % Thrust

                drag(i,1)               = 0.5 * rho_air(i,1) * (vel((i-1),1)).^2 * cd_curve(Ma(i-1,1)) * Ar;                %Total drag

                Ay(i,1)                 = (thrust(i-1,1)*sind(obj.input.sim.thetaL) - drag((i-1),1)*sind(obj.input.sim.thetaL) + mtot(i-1,1)*g)./mtot(i-1,1);   % Acceleration in y direction(m/s)

                Ax(i,1)                 = (thrust(i-1,1)*cosd(obj.input.sim.thetaL) - drag((i-1),1)*cosd(obj.input.sim.thetaL))./mtot(i-1,1);% Acceleration in x direction(m/s)

                gforce(i,1)             = sqrt(Ay(i,1)^2 + Ax(i,1)^2)/-g;                                                   % Acceleration (G)

                force(i,1)              = mtot(i, 1)*(sqrt(Ay(i,1)^2 + Ax(i,1)^2)) + drag(i,1);                             % Total applied axial load

                velx(i,1)              = velx((i-1),1) + Ax(i,1)*obj.input.sim.dt;                                          % Velocity in x direction at time t

                vely(i,1)              = vely((i-1),1) + Ay(i,1)*obj.input.sim.dt;                                          % Velocity in y direction at time t

                vel(i,1)                = sqrt(velx(i,1).^2 + vely(i,1).^2);                                                % Average velocity at time t

                sy(i,1)                 = sy((i-1),1) + vely((i-1),1) * obj.input.sim.dt + (0.5*Ay(i-1,1)*obj.input.sim.dt.^2); % Height at time t

                sx(i,1)                 = sx((i-1),1) + velx((i-1),1)*obj.input.sim.dt + (0.5*Ax(i-1,1)*obj.input.sim.dt.^2);   % x distance at time t

                Ma(i,1)                 = vel(i,1)/(sqrt(1.4*R_air*T_curve(sy(i-1,1))));                                        % Mach number at time t

            end

            prop.performance.thrust = thrust;
            prop.performance.isp    = Isp;
            prop.performance.Pcc    = Pcc;
            prop.performance.Tcc    = Tcc;

            obj.flight.y            = sy;
            obj.flight.x            = sx;
            obj.flight.vely         = vely;
            obj.flight.velx         = velx;
            obj.flight.Ma           = Ma;
            obj.flight.gforce       = gforce;
            obj.flight.force        = force;
            obj.flight.mtot         = mtot;
            obj.flight.drag         = drag;

        end

        function obj    = getChamberPressure(obj)
        end

        function obj    = getChamberTemperature(obj)
        end

        function obj    = simulate(obj)
        end

    end

end
