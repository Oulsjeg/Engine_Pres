function [dn,m_dot, dT, P] = C2oolPropTest(n, T, V, m_tank, A_inj, Cd, Pe, nPres)
%OX_TANK_SIM Models thermochemical processes of blow-down in oxidizer tank.
%   INPUTS:
%       "n" = 1x2 vector of ox moles in tank (kmol).
%       "T" = ox tank temperature (K).
%       "V" = volume of ox tank (m^3).
%       "m_tank" = mass of oxidizer tank (kg).
  %    "m_dot_prev" = ox mass flow rate of previous time step (kg/s).
%       "A_inj" = total orifice area of injector plate (m^2).
%       "Cd" = injector plate discharge coefficient.
%       "Pe" = pressure at exit of injector plate (Pa).
%   OUTPUTS:
%       "m_dot" = scalar ox mass flow rate (kg/s).
%       "dn" = 1x2 differential change in ox moles (kmol).
%       "dT" = differential change in ox tank temperature (K).
%       "P" = scalar internal pressure of ox tank (Pa).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V=0.03;
n=[250 125];
m_tank=10;
A_inj=0.02891;
Cd=0.65;
Pe= 20000;
nPres=50;

R  = 8314;
MW = 44.013;   %Molar mass of N2O (kg/kmol).
T=298;
            

                % Molar specific vol. of liq. N2O [m^3/kmol]
                Vhat_l      = MW / py.CoolProp.CoolProp.PropsSI('D', 'P', 100000, 'T', T(1), 'N2O');

                % Molar c_V of Pressurant gas [J/(kmol*K)]
                CVhat_Pres  =MW * py.CoolProp.CoolProp.PropsSI('O', 'P', 100000, 'T', T(1));

                % Molar c_V of N2O gas [J/(kmol*K)]
                CVhat_g     = MW * py.CoolProp.CoolProp.PropsSI('O', 'T', T(1), 'Q', 1, 'N2O') - R;

                % Molar c_V of N2O liq approx. = c_P [J/(kmol*K)]
                CVhat_l     = MW * py.CoolProp.CoolProp.PropsSI('C', 'T', T(1), 'Q', 1, 'N2O');

                % Heat of vaporization of N2O [J/kmol]
                Hv          = py.CoolProp.CoolProp.PropsSI('H', 'T', T(1), 'Q', 1, 'N2O');
                Hl          = py.CoolProp.CoolProp.PropsSI('H', 'T', T(1), 'Q', 0, 'N2O');
                delta_Hv    = MW * (Hv - Hl);

                % Vapour P of N20 (Pa).
                P_sat       = py.CoolProp.CoolProp.PropsSI('P', 'T', T(1), 'Q', 0, 'N2O');

                % Derivative of vapour P with respect to T.
                dP_sat      = py.CoolProp.CoolProp.PropsSI('d(P)/d(T)|D', 'T', T, 'Q', 1, 'N2O');
               
                c_P         = (4.8 + 0.00322*T)*155.239;                    %Specific heat of tank, Aluminum [J/kg*K]


            P           = (nPres + n(2))*R*T/(V - n(1)*Vhat_l);
            a           = m_tank*c_P + nPres*CVhat_Pres + n(2)*CVhat_g + n(1)*CVhat_l;
            RT          = P*Vhat_l;
            e           = -delta_Hv + RT;
            f           = -Cd*A_inj*sqrt(abs(2/MW*(P - Pe)/Vhat_l));
            j           = -Vhat_l*P_sat;
            k           = (V- n(1)*Vhat_l)*dP_sat;
            l           = R*T;
            q           = R*n(2);

            Z           = (-f*(-j*a + (q - k)*RT))/(a*(l + j) + (q - k)*(e - RT));
            W           = (-Z*(l*a + (q - k)*e))/(-j*a + (q - k)*RT);

            dT          = (RT*W + e*Z)/a;
            


%Time derivative of tank conditions
            m_dot     = -f*MW;         % Mass flow rate
            dn       = [W, Z];                       % Time derivative of nitrous amount
            dT       = [dT, dT];                     % Time derivative of tank temperature

end

