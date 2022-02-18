clear;

V       = 0.0025;         % Volume    [m^3]
Ti      = 310;          % Temp      [K]
cnv     = 6894.76;
Pi      = 2500*cnv;     % Pressure  [psi]
gas     = 'N2';

rho     = py.CoolProp.CoolProp.PropsSI('D', 'T', Ti, 'P', Pi, gas);

time    = 15;
num     = 300;
t       = linspace(0, time, num)';
dt      = time/(num+1);

P       = zeros(num, 1);
R       = zeros(num, 1);
T       = zeros(num, 1);
m       = zeros(num, 1);
dPdt    = zeros(num, 1);
dTdt    = zeros(num, 1);

mdot    = linspace(0, 0.01, num)';
R(1, 1) = rho;
T(1, 1) = Ti;
m(1, 1) = V*rho;
P(1, 1) = Pi;
qdot    = linspace(0, 1000, num)';

ru      = 8314 / 28;

for i = 2:num

    [dPdt_i, dTdt_i, dZdt]  = n2_blowdown(V, T(i-1, 1), P(i-1, 1), mdot(i-1, 1), qdot(i-1, 1), m(i-1, 1), gas);
    dPdt(i, 1)      = dPdt_i;
    dTdt(i, 1)      = dTdt_i;

    P(i, 1)     = P(i-1, 1) + dPdt(i, 1)*dt;
    T(i, 1)     = T(i-1, 1) + dTdt(i, 1)*dt;
    m(i, 1)     = m(i-1, 1) - mdot(i, 1)*dt;
    rho         = P(i, 1) / (ru * getZ(T(i, 1), P(i, 1)/cnv) * T(i, 1));
    R(i, 1)     = py.CoolProp.CoolProp.PropsSI('D', 'P', P(i, 1), 'T', T(i, 1), gas);

end

P = P/cnv;

figure(1)
plot(t, P, 'Color', 'black');
title('Pressure (psi)');

figure(2)
plot(t, T, 'Color', 'black');
title('Temperature [K]');


