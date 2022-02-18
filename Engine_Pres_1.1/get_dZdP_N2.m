function dZdP = get_dZdP_N2(T, P)
    
    P           = P / 6894.76;
    
    nitrogen.h.theta    = [0.90370032155133;...
                          -3.99164830787538;...
                           4.44554878990612;...
                          -1.68387394930366;...
                           1.84282855081908;...
                          -2.71807522455834;...
                           1.80658523674363;...
                          -0.00026662830718;...
                           0.16405364316350];

    nitrogen.h.alpha    = [1;               ...
                           0.45607085009281;...
                           0.99224794564113;...
                           1.58495789262624;...
                           0.53133147588636;...
                           1.29132167947510;...
                           1.44008913900161;...
                           2.74997487910292;...
                           2.36611999082672];

    nitrogen.l.theta    = [0.46742656471647;...
                          -0.53799565472298;...
                          -9.22454428760102;...
                           9.15603503101003;...
                           3.18808664459882;...
                           0.30163700042055;...
                          -0.27300234680706;...
                          -1.00749719408221;...
                          -1.49106816983329];

    nitrogen.l.alpha    = [1;               ...
                           1.41102397459172;...
                           0.33562799290636;...
                           0.79810083070486;...
                           0.01008992455881;...
                           2.53968667359886;...
                           2.51281397715323;...
                           1.20879498088509;...
                           1.69572064361084];
    
    if (T >= 220 && T <= 320) && (P >= 800 && P <= 4500)

        nh          = nitrogen.h;
    	T           = T / 126.2;
        P           = P / 492.3;

        dZdPr       = 0.7*nh.theta(2)*T^nh.alpha(2)*P^(-0.3)        ...
                    + 0.5*nh.theta(3)*T^nh.alpha(3)*P^(-0.5)        ...
                    + 0.3*nh.theta(4)*T^nh.alpha(4)*P^(-0.7)        ...
                    + nh.theta(5)*T^nh.alpha(5)                     ...
                    + nh.theta(6)*T^nh.alpha(6)                     ...
                    + nh.theta(7)*T^nh.alpha(7)                     ...
                    + nh.theta(8)*nh.alpha(8)*P^(nh.alpha(8) - 1);

        dZdP        = dZdPr / 492.3;

    elseif (T >= 180 && T < 220) && (P >= 800 && P <= 4500)

        nl          = nitrogen.l;
    	T           = T / 126.2;
        P           = P / 492.3;

        dZdPr       = 0.7*nl.theta(2)*T^nl.alpha(2)*P^(-0.3)        ...
                    + 0.5*nl.theta(3)*T^nl.alpha(3)*P^(-0.5)        ...
                    + 0.3*nl.theta(4)*T^nl.alpha(4)*P^(-0.7)        ...
                    + nl.theta(5)*T^nl.alpha(5)                     ...
                    + nl.theta(6)*T^nl.alpha(6)                     ...
                    + nl.theta(7)*T^nl.alpha(7)                     ...
                    + nl.theta(8)*nl.alpha(8)*P^(nl.alpha(8) - 1);

        dZdP        = dZdPr / 492.3;

    else
        if (T < 180 || T > 320)
            error('utilitiesClass > get_dZdP(): Temperature is out off acceptable bounds (T = %f)', T);
        elseif (P < 800 || P > 4500)
            error('utilitiesClass > get_dZdP(): Pressure is out off acceptable bounds (P = %f)', P);
        end
    end