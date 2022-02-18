classdef utilitiesClass < handle

    properties (Constant)
        g0          = 9.8056;
        earthR      = 6371e3;
        airR        = 287.05;
        R           = 8314;
        airM        = 0.0289644; %molar mass of air (kg/mol)
    end

    properties
        zeroArray
        noxProp
        nitrogen
        T_curve
    end

    methods    (Static)
        function obj                = utilitiesClass(input)

            obj.zeroArray           = zeros(input.sim.numpt, 1);

            obj.noxProp.rho_c       = 452;
            obj.noxProp.Tc          = 309.57;
            obj.noxProp.Pc          = 7.251e6;

            obj.noxProp.Coefs.V1    = 96.512;
            obj.noxProp.Coefs.V2    = -4045;
            obj.noxProp.Coefs.V3    = -12.277;
            obj.noxProp.Coefs.V4    = 2.886e-5;
            obj.noxProp.Coefs.V5    = 2;

            obj.noxProp.Coefs.T1    = 2.3215e7;
            obj.noxProp.Coefs.T2    = 0.384;
            obj.noxProp.Coefs.T3    = 0;
            obj.noxProp.Coefs.T4    = 0;

            obj.noxProp.Coefs.Q1    = 2.781;
            obj.noxProp.Coefs.Q2    = 0.27244;
            obj.noxProp.Coefs.Q3    = 309.57;
            obj.noxProp.Coefs.Q4    = 0.2882;

            obj.noxProp.Coefs.D1    = 0.2934e5;
            obj.noxProp.Coefs.D2    = 0.3236e5;
            obj.noxProp.Coefs.D3    = 1.1238e3;
            obj.noxProp.Coefs.D4    = 0.2177e5;
            obj.noxProp.Coefs.D5    = 479.4;

            obj.noxProp.Coefs.E1    = 6.7556e4;
            obj.noxProp.Coefs.E2    = 5.4373e1;
            obj.noxProp.Coefs.E3    = 0;
            obj.noxProp.Coefs.E4    = 0;
            obj.noxProp.Coefs.E5    = 0;

            obj.noxProp.Coefs.b1    = 1.72328;
            obj.noxProp.Coefs.b2    = -0.8395;
            obj.noxProp.Coefs.b3    = 0.5106;
            obj.noxProp.Coefs.b4    = 0.10412;

            obj.noxProp.Coefs.q1    = -6.71893;
            obj.noxProp.Coefs.q2    = 1.35966;
            obj.noxProp.Coefs.q3    = 1.3779;
            obj.noxProp.Coefs.q4    = -4.051;

            obj.nitrogen.Coefs.C1    = 0.28883e5;
            obj.nitrogen.Coefs.C2    = 0;
            obj.nitrogen.Coefs.C3    = 0;
            obj.nitrogen.Coefs.C4    = 0;
            obj.nitrogen.Coefs.C5    = 0;

        end
    end

    methods

        function output = getN2OProperties(obj, T)
            
            T = T/obj.noxProp.Tc;
            
            A   = obj.noxProp.Coefs.b1*(1-T)^(1/3)  + ...
                  obj.noxProp.Coefs.b2*(1-T)^(2/3)  + ...
                  obj.noxProp.Coefs.b3*(1-T)        + ...
                  obj.noxProp.Coefs.b4*(1-T)^(4/3);

            Ap  = obj.noxProp.Coefs.q1*(1-T)        + ...
                  obj.noxProp.Coefs.q2*(1-T)^(3/2)  + ...
                  obj.noxProp.Coefs.q3*(1-T)^(5/2)  + ...
                  obj.noxProp.Coefs.q4*(1-T)^(5);

            output.rho  = obj.noxProp.rho_c*exp(A);
            output.P    = obj.noxProp.Pc*exp(Ap/T);
        end

        function rho = getRho(~, T, P, name)
            
            if strcmp(name, 'C2H5OH')
                name = 'ethanol';
            end
            
            rho = py.CoolProp.CoolProp.PropsSI('D', 'P', P, 'T', T, name);
        end

        function Pamb = getAmbientPressure(~, altitude)

            Ru = obj.R/1000;%universal gas constant (J/mol/K)

            if ( 0 <= altitude) && (altitude < 11000)
                Lb = -0.0065; %Lapse rate (K/m)
                Tb = 288.15; %Standard temperature (K)
                Pb = 101325; %Static pressure (Pa)
                h0 = 0;
                Pamb = Pb * (Tb/(Tb + Lb*(altitude - h0))) ^ (obj.g0*obj.M/Ru/Lb);

            elseif ( 11000 <= altitude) && (altitude < 20000)
                Tb = 216.65;
                Pb = 22632.1;
                h1 = 11000;
                Pamb = Pb * exp(-obj.g0*obj.M*(altitude - h1)/Ru/Tb);

            elseif ( 20000 <= altitude) && (altitude < 32000)
                Lb = 0.001;
                Tb = 216.65;
                Pb = 5474.89;
                h2 = 20000;
                Pamb = Pb * (Tb/(Tb + Lb*(altitude - h2))) ^ (obj.g0*obj.M/Ru/Lb);

            elseif (32000 <= altitude) && (altitude < 47000)
                Lb = 0.0028;
                Tb = 228.65;
                Pb = 868.02;
                h3 = 32000;
                Pamb = Pb * (Tb/(Tb + Lb*(altitude - h3))) ^ (g*obj.M/Ru/Lb);


            elseif(47000 <= altitude) && (altitude < 51000)
                Tb = 270.65;
                Pb = 110.91;
                h4 = 47000;
                Pamb = Pb * exp(-obj.g0*obj.M*(altitude - h4)/Ru/Tb);

            elseif(51000 <= altitude) && (altitude < 71000)
                Lb = -0.0028;
                Tb = 270.65;
                Pb = 66.94;
                h5 = 51000;
                Pamb = Pb * (Tb/(Tb + Lb*(altitude - h5))) ^ (obj.g0*obj.M/Ru/Lb);

            elseif(71000 <= altitude) && (altitude < 86000)
                Lb = -0.002;
                Tb = 214.65;
                Pb = 3.96;
                h6 = 71000;
                Pamb = Pb * (Tb/(Tb + Lb*(altitude - h6))) ^ (obj.g0*obj.M/Ru/Lb);

            else
                Pamb = 0;
            end
        end
        
        function getTcurve(obj)
            
            altitude                = [0,11000,20000,32000,47000,51000,71000];
            temps                   = [288,216,216,228,270,270,214];
            obj.T_curve             = fit(altitude',temps','cubicinterp');
        end
    end
end
