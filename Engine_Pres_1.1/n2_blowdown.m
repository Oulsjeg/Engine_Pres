function [dPdt, dTdt, dZdt] = n2_blowdown(V, T, P, mdot, qdot, m, gas)

    R = 8314 / 28;
    cnv = 6894.76;

    Zt = get_dZdT_N2(T, P);
    Zp = get_dZdP_N2(T, P)/cnv;
    Z = py.CoolProp.CoolProp.PropsSI('Z', 'P', P, 'T', T, gas);
    cv = py.CoolProp.CoolProp.PropsSI('O', 'P', P, 'T', T, gas);

    dPdt = (m*R*Z^2 + P*V*Zt)*qdot/(V*cv*(m*(Z - P*Zp)))    ...
         - mdot*((cv + Z*R)*m*Z + P*V*Zt)*(P/(m*cv))/(m*(Z - P*Zp));

    dZdt = (m*R*Z^2 * Zp + V*Z*Zt)*dPdt/(m*R*Z^2 + P*V)     ...
         - Zt*(P*V*Z/m)*mdot/(m*R*Z^2 + P*V);

    dTdt = V*(dPdt/(m*Z) - dZdt*P/(m*Z^2) - mdot*P/(m*Z^2))/R;
end