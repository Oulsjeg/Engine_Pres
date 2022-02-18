classdef propellantTankClass < handle

    properties (Access = public)

        name            % Propellant name                           (str)
        type            % Propellant type                           (str)
        inputTag        % Tag for name of dedicated input field     (str)
        input           % Input structure                           (struct)
        tank            % Tank body                                 (struct)
        propellant      % Liquid phase                              (fluidClass)
        pressurant      % Pressurant vapor phase                    (fluidClass)
        cg              % CG of propellant tank and contents        (double)
        m               % Mass of propellant tank and contents      (double)
        ullage          % Ullage volume                             (double)
        offset          % Distance offset before tank system        (double)

        isPressurized   % Pressurization flag                       (bool)
        pressurantOrder % Pressurant tank arrangment                (bool)
        utilities       % Utilities                                 (utilitiesClass)
    end

    methods

        %-----------------------------------------------------------------------
        %   METHOD: propellantTankClass
        %   Constructs propellantTankClass given required inputs.
        %
        %   INPUTS .............................................................
        %     - <input> (struct): A structure that contains all input variables
        %     - <inputTag> (str): A tag that identifies the field name in
        %                         <input> that contains data relevant to the
        %                         specifid tank being created.
        %   OUTPUTS ............................................................
        %     - <obj> (class):    Returns created propellantTankClass
        %-----------------------------------------------------------------------
        function obj    = propellantTankClass(inputTag, input)

            % Populate general properties
            obj.input           = input;                            % Input
            obj.inputTag        = inputTag;
            obj.name            = input.(inputTag).name;            % Propellant name
            obj.type            = input.(inputTag).fluidtype;       % Propellant type
            obj.m               = zeros(input.sim.numpt, 1);        % Mass
            obj.cg              = zeros(input.sim.numpt, 1);        % CG

            % Create utilities property
            obj.utilities       = utilitiesClass(input);

            % Populate tank struct
            obj.tank.t          = input.design.tankThickness;       % Wall thickness
            obj.tank.r          = input.design.diameter/2;          % Radius
            obj.tank.m          = input.(inputTag).mTank;           % Mass
            obj.tank.l          = input.(inputTag).lTank;           % length
            obj.tank.v          = input.(inputTag).vTank;           % Volume
            obj.tank.offset     = input.(inputTag).offset;          % Offset distance
            obj.tank.cg         = obj.tank.l/2;                     % CG location

            initP               = input.(inputTag).Pinit;
            initT               = input.(inputTag).Tinit;
            initm               = input.(inputTag).mInit;

            obj.tank.ullage     = obj.tank.v                                    ...
                                - initm/obj.utilities.getRho(initT, initP, input.(inputTag).name);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Create propellant liquid, vapor and pressurant vapor classes

            % Initialize propellant liquid phase
            obj.propellant      = twoPhaseFluidClass(input, inputTag);      % Liquid phase

            initStruct.m        = [input.(inputTag).mInit, 0];
            initStruct.T        = [input.(inputTag).Tinit, input.(inputTag).Tinit];
            initStruct.P        = [input.(inputTag).Pinit, input.(inputTag).Pinit];

            obj.propellant.setInitialConditions(initStruct);

            % Initialize pressurant vapor phase
            obj.pressurant      = onePhaseFluidClass(input, input.(inputTag).pressurant);      % Pres vapor phase

            initStruct.T        = input.(inputTag).Tinit;
            initStruct.P        = input.(inputTag).Pinit;
            initStruct.m        = obj.tank.ullage                           ...
                                * obj.utilities.getRho(initStruct.T, initStruct.P, obj.pressurant.name);

            obj.pressurant.setInitialConditions(initStruct);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.isPressurized   = input.(inputTag).isPressurized;   % isPressurized flag

            % Select blowdown function
            obj.setBlowdownCharacteristics();
        end

        %-----------------------------------------------------------------------
        %   METHOD: simulateBlowdown
        %   Selects and calls appropriate blowdown function depending on the
        %   type of blowdown specified in the input file.
        %
        %   INPUTS: NONE
        %   OUTPUTS: NONE
        %-----------------------------------------------------------------------
        function          getBlowdown(obj)

            switch obj.type
                case 'Oxidizer'
                    switch obj.input.(obj.inputTag).blowdownMode
                        case 'constantPressure'
                            obj.OX_constantPressureBlowdown();
                        otherwise
                            error('Blowdown mode %s has not been specialized', obj.input.design.presMode);
                    end
                    
                case 'Fuel'
                    switch obj.input.(obj.inputTag).blowdownMode
                        case 'constantMdot'
                            obj.FUEL_constantMdotBlowdown();
                        otherwise
                            error('Blowdown mode %s has not been specialized', obj.input.design.presMode);
                    end   
                    
                otherwise
                    error('Propellant type %s unknown. Type can either be Oxidizer or Fuel', obj.type);
            end
        end

        %-----------------------------------------------------------------------
        %   METHOD: setBlowdownFunction
        %   Selects the functions that computes the blowdown characteristics for
        %   the given propellant.
        %
        %   INPUTS: NONE
        %   OUTPUTS: NONE
        %-----------------------------------------------------------------------
        function          setBlowdownCharacteristics(obj)

            switch obj.type
                case 'Oxidizer'
                    if strcmp(obj.name, 'N2O')
                        obj.propellant.blowdownCharacteristics  = @ obj.bdChars_NIROUSOXIDE;
                    else
                        error('Oxidizer type %s has no mass flow rate function specified', obj.name);
                    end
                    
                case 'Fuel'
                    obj.propellant.blowdownCharacteristics      = @ obj.bdChars_LIQUIDS_MDOT;
                        
                otherwise
                    error('Blowdown function for propellant type %s has not been specialized.', obj.type);
            end
        end

        %-----------------------------------------------------------------------
        %   METHOD: constantPressureBlowdown
        %   Performs propellant blowdown using the constant tank pressure
        %   pressurization mode. Uses the selected blowdown function from
        %   setBlowdownFunction() to compute the blowdown characteristics at
        %   each time step.
        %
        %   INPUTS: NONE
        %   OUTPUTS: NONE
        %-----------------------------------------------------------------------
        function          OX_constantPressureBlowdown(obj)

            in.mDot             = obj.input.design.mDotox;
            in.T                = obj.input.(obj.inputTag).Tinit;
            
            in.V                = obj.tank.v;
            in.mTank            = obj.tank.m;
            in.P                = obj.propellant.IC.P;
            in.n                = obj.propellant.IC.n;
            in.Cd               = obj.input.design.injCd;
            in.Pe               = obj.input.design.Pcc  * obj.input.settings.cnv;
            in.nPres            = obj.pressurant.m(1, 1)/obj.pressurant.MW;
            in.model            = obj.input.settings.fluidModel;

            in                  = obj.findAinj(in);

            for i = 2 : obj.input.sim.numpt

                out             = obj.propellant.blowdownCharacteristics(in);
                dP              = obj.propellant.P(i - 1, 1) - out.P;

                while (dP > obj.input.settings.dPtol)

                    in.nPres        = in.nPres  + obj.input.settings.npr_inc;
                    out             = obj.propellant.blowdownCharacteristics(in);
                    dP              = obj.propellant.P(i - 1, 1)   - out.P;
                end

                obj.propellant.n(i, :)      = obj.propellant.n(i-1, :) + out.dn * obj.input.sim.dt;
                obj.propellant.T(i, :)      = obj.propellant.T(i-1, 1) + out.dT * obj.input.sim.dt;
                obj.propellant.mdot(i, 1)   = out.mDot;
                obj.propellant.P(i, :)      = [out.P, out.P];
                
                obj.pressurant.n(i, 1)      = in.nPres;

                in.T                        = obj.propellant.T(i, 1);
                in.n                        = obj.propellant.n(i, :);
                in.P                        = obj.propellant.P(i, 1);
            end

            obj.propellant.m        = obj.propellant.n * obj.propellant.MW;
            obj.pressurant.m        = obj.pressurant.n * obj.pressurant.MW;
        end
        
        function          FUEL_constantMdotBlowdown(obj)
            
            mdot    = obj.input.design.mDotox / obj.input.design.OF;
            T       = obj.input.(obj.inputTag).Tinit;
            
            obj.propellant.mdot(1, 1) = mdot;
            
            for i = 2 : obj.input.sim.numpt
                
                out = obj.bdChars_LIQUIDS_MDOT(mdot);
                obj.propellant.m(i, 1)      = obj.propellant.m(i-1, 1) + out.mDot * obj.input.sim.dt;
                obj.propellant.mdot(i, 1)   = mdot;
                obj.propellant.T(i, 1)      = T;
            end
            
        end
        
        function output = bdChars_LIQUIDS_MDOT(~, mdot)
            output.mDot = -mdot;
        end

        % Nitrous oxide blowdown
        function output = bdChars_NIROUSOXIDE(obj, input)

            % Input: T, nPres, nVap, nLiq, V, mTank, Cd, Ainj, Pe
            if strcmp(input.model, 'ideal')

                % Molar specific vol. of liq. N2O [m^3/kmol]
                Vhat_l      = obj.utilities.noxProp.Coefs.Q2                          ^ ...
                             (1 + (1 - input.T/obj.utilities.noxProp.Coefs.Q3)        ^ ...
                              obj.utilities.noxProp.Coefs.Q4)                         / ...
                              obj.utilities.noxProp.Coefs.Q1;

                % Molar c_V of He [J/(kmol*K)]
                CVhat_Pres  = obj.utilities.nitrogen.Coefs.C1                   + ...
                              obj.utilities.nitrogen.Coefs.C2*input.T           + ...
                              obj.utilities.nitrogen.Coefs.C3*input.T^2         + ...
                              obj.utilities.nitrogen.Coefs.C4*input.T^3         + ...
                              obj.utilities.nitrogen.Coefs.C5*input.T^4         - obj.utilities.R;

                a           = obj.utilities.noxProp.Coefs.D3/input.T;
                b           = obj.utilities.noxProp.Coefs.D5/input.T;

                % Molar c_V of N2O gas [J/(kmol*K)]
                CVhat_g     = obj.utilities.noxProp.Coefs.D1                    + ...
                              obj.utilities.noxProp.Coefs.D2*(a/sinh(a))^2      + ...
                              obj.utilities.noxProp.Coefs.D4*(b/cosh(b))^2      - obj.utilities.R;

                % Molar c_V of N2O liq approx. = c_P [J/(kmol*K)]
                CVhat_l     = obj.utilities.noxProp.Coefs.E1                    + ...
                              obj.utilities.noxProp.Coefs.E2*input.T            + ...
                              obj.utilities.noxProp.Coefs.E3*input.T^2          + ...
                              obj.utilities.noxProp.Coefs.E4*input.T^3          + ...
                              obj.utilities.noxProp.Coefs.E5*input.T^4;

                % Reduced temperature.
                Tr          = input.T/obj.utilities.noxProp.Tc;

                % Heat of vaporization of N2O [J/kmol]
                delta_Hv    = obj.utilities.noxProp.Coefs.T1*(1 - Tr)           ^ ...
                             (obj.utilities.noxProp.Coefs.T2                    + ...
                              obj.utilities.noxProp.Coefs.T3*Tr                 + ...
                              obj.utilities.noxProp.Coefs.T4*Tr^2);

                % Vapour P of N20 (Pa).
                P_sat       = exp(obj.utilities.noxProp.Coefs.V1                + ...
                              obj.utilities.noxProp.Coefs.V2/input.T            + ...
                              obj.utilities.noxProp.Coefs.V3*log(input.T)       + ...
                              obj.utilities.noxProp.Coefs.V4*input.T^obj.utilities.noxProp.Coefs.V5);

                % Derivative of vapour P with respect to T.
                dP_sat      = (-obj.utilities.noxProp.Coefs.V2/(input.T^2)                      + ...
                              obj.utilities.noxProp.Coefs.V3/input.T                            + ...
                              obj.utilities.noxProp.Coefs.V4*obj.utilities.noxProp.Coefs.V5     * ...
                              input.T^(obj.utilities.noxProp.Coefs.V5-1))                       * ...
                              exp(obj.utilities.noxProp.Coefs.V1                                + ...
                              obj.utilities.noxProp.Coefs.V2/input.T                            + ...
                              obj.utilities.noxProp.Coefs.V3*log(input.T)                       + ...
                              obj.utilities.noxProp.Coefs.V4*input.T^obj.utilities.noxProp.Coefs.V5);

            elseif strcmp(input.model, 'coolprop')
                % Molar specific vol. of liq. N2O [m^3/kmol]
                Vhat_l      = obj.propellant.MW / py.CoolProp.CoolProp.PropsSI('D', 'P', input.P(1), 'T', input.T(1), 'N2O');

                % Molar c_V of Pressurant gas [J/(kmol*K)]
                CVhat_Pres  = obj.pressurant.MW * py.CoolProp.CoolProp.PropsSI('O', 'P', input.P(1), 'T', input.T(1), obj.pressurant.name);

                % Molar c_V of N2O gas [J/(kmol*K)]
                CVhat_g     = obj.propellant.MW * py.CoolProp.CoolProp.PropsSI('O', 'T', input.T(1), 'Q', 1, 'N2O') - obj.utilities.R;

                % Molar c_V of N2O liq approx. = c_P [J/(kmol*K)]
                CVhat_l     = obj.propellant.MW * py.CoolProp.CoolProp.PropsSI('C', 'T', input.T(1), 'Q', 1, 'N2O');

                % Heat of vaporization of N2O [J/kmol]
                Hv          = py.CoolProp.CoolProp.PropsSI('H', 'T', input.T(1), 'Q', 1, 'N2O');
                Hl          = py.CoolProp.CoolProp.PropsSI('H', 'T', input.T(1), 'Q', 0, 'N2O');
                delta_Hv    = obj.propellant.MW * (Hv - Hl);

                % Vapour P of N20 (Pa).
                P_sat       = py.CoolProp.CoolProp.PropsSI('P', 'T', input.T(1), 'Q', 0, 'N2O');

                % Derivative of vapour P with respect to T.
                dP_sat      = py.CoolProp.CoolProp.PropsSI('d(P)/d(T)|D', 'T', input.T, 'Q', 1, 'N2O');
            else
                error('wat')
            end

            %Specific heat of tank, Aluminum [J/kg*K]
            c_P         = (4.8 + 0.00322*input.T)*155.239;

            %Simplified expression definitions for solution:
            P           = (input.nPres + input.n(2))*obj.utilities.R*input.T/(input.V - input.n(1)*Vhat_l);
            a           = input.mTank*c_P + input.nPres*CVhat_Pres + input.n(2)*CVhat_g + input.n(1)*CVhat_l;
            RT          = P*Vhat_l;
            e           = -delta_Hv + obj.utilities.R*input.T;
            f           = -input.Cd*input.Ainj*sqrt(abs(2/obj.propellant.MW*(P - input.Pe)/Vhat_l));
            j           = -Vhat_l*P_sat;
            k           = (input.V- input.n(1)*Vhat_l)*dP_sat;
            l           = obj.utilities.R*input.T;
            q           = obj.utilities.R*input.n(2);

            Z           = (-f*(-j*a + (q - k)*RT))/(a*(l + j) + (q - k)*(e - RT));
            W           = (-Z*(l*a + (q - k)*e))/(-j*a + (q - k)*RT);

            dT          = (RT*W + e*Z)/a;

            %Time derivative of tank conditions
            output.mDot     = -f*obj.propellant.MW;         % Mass flow rate
            output.dn       = [W, Z];                       % Time derivative of nitrous amount
            output.dT       = [dT, dT];                     % Time derivative of tank temperature
            output.P        = P;
        end        
        
        % No sure if this is still needed
        function          setupTankSystem(obj, data)

            obj.liq.m           = data.mass.l;
            obj.propVap.m       = data.mass.vo;
            obj.presVap.m       = data.mass.vp;

            obj.liq.mdot        = data.mdot.l;
            obj.propVap.mdot    = data.mdot.vo;
            obj.presVap.mdot    = data.mdot.vp;

            obj.liq.rho         = data.rho;

            obj.mprop           = obj.tank.m.*ones(length(obj.liq.m), 1)    ...
                                + obj.liq.m                                 ...
                                + obj.propVap.m                             ...
                                + obj.presVap.m;

            if obj.isPressurized
                obj.pressurant.setupTankSystem(obj.presVap.mdot);
            end

        end

        % Get tank system CG
        function          getCG(obj)

            obj.liq.l(1)         = obj.liq.m(1)/(obj.liq.rho(1)*pi*(obj.tank.r - obj.tank.t)^2);
            obj.liq.cg(1)        = obj.liq.l(1)/2;

            obj.propVap.l(1)      = obj.tank.l                - obj.liq.l(1);
            obj.propVap.cg(1)     = obj.propVap.l(1)/2;

            obj.presVap.l(1)      = obj.tank.l                - obj.liq.l(1);
            obj.presVap.cg(1)     = obj.presVap.l(1)/2;

            obj.cgprop(1)           = (obj.liq.m(1)             *   (obj.propVap.l(1) + obj.liq.cg(1))  ...
                                    +  obj.propVap.m(1)         *   obj.presVap.cg(1)                   ...
                                    +  obj.presVap.m(1)         *   obj.presVap.cg(1)                   ...
                                    +  obj.tank.m*obj.tank.cg)  /   obj.mprop(1);

            for i   = 2:length(obj.mprop)

                dl                  = (obj.liq.m(i) - obj.liq.m(i - 1)) ...
                                    / (pi*obj.liq.rho(i - 1)*(obj.tank.r - obj.tank.t)^2);

                obj.liq.cg(i)       = obj.liq.cg(i - 1)         + dl/2;
                obj.liq.l(i)        = obj.liq.l(i - 1)          + dl;

                obj.propVap.cg(i)   = obj.propVap.cg(i - 1)     - dl/2;
                obj.propVap.l(i)    = obj.propVap.l(i - 1)      - dl;

                obj.presVap.cg(i)   = obj.presVap.cg(i - 1)     - dl/2;
                obj.presVap.l(i)    = obj.presVap.l(i - 1)      - dl;

                obj.cgprop(i)       = (obj.liq.m(i)             *   (obj.propVap.l(i) + obj.liq.cg(i))  ...
                                    +  obj.propVap.m(i)         *   obj.presVap.cg(i)                   ...
                                    +  obj.presVap.m(i)         *   obj.presVap.cg(i)                   ...
                                    +  obj.tank.m*obj.tank.cg)  /   obj.mprop(i);
            end
        end

        % Find required injector orifice area to supply desired mass flow rate
        function input  = findAinj(obj, input)

            out.mDot    = 1000;
            input.Ainj  = 0.0001;

            while (out.mDot > input.mDot)

                input.Ainj          = input.Ainj - 0.0000001;
                out                 = obj.propellant.blowdownCharacteristics(input);
            end

            input.mDot = out.mDot;
        end
    end

end
