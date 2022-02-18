%--------------------------------------------------------------------------
% AUTHORS:
%   1) Muhammad Bonse, Ryerson University
%--------------------------------------------------------------------------

function [engine] = combustion()

% Define Problem Parameters
OF_i        = 1.5;  % Initial OF
OF_f        = 3;    % Final OF
nValues     = 10;   % Number of data points to generate

OF      =   linspace(OF_i,OF_f,nValues);      % OF Ratios
num_OF  =   length(OF);                       % Number of OF ratios
Pcc     =   100000;                       % Chamber Pressure [Psi]
PcPe    =  100; % Chamber/Exit Pressure. Get from Nozzle Code

fuel    =   {'C2H5OH', 100, 298};
ox      =   {'N2O', 100, 278};
    
%[OF, num_OF, Pcc, PcPe, fuel, ox] = tables(OF_i, OF_f, nValues, rocket.Pcc_des, rocket.alt_Bout);

% Loop through all OF ratios and solve
[CEA_result, cstar,Tcc,Te,Pe,gamm_e,Isp,Ma_e,rho_e,R_e] = ...
    CEA_tables_solve(num_OF, OF, Pcc, 10, PcPe, fuel, ox);

% Generating curves for cea results
engine.cstar_curve = fit(OF',cstar','cubicinterp');
engine.Tcc_curve = fit(OF',Tcc','cubicinterp');
engine.Te_curve = fit(OF',Te','cubicinterp');
engine.Pe_curve = fit(OF',Pe','cubicinterp');
engine.gamm_e_curve = fit(OF',gamm_e','cubicinterp');
engine.Isp_curve = fit(OF',Isp','cubicinterp');
engine.Ma_e_curve = fit(OF',Ma_e','cubicinterp');
engine.rho_e_curve = fit(OF',rho_e','cubicinterp');
engine.R_e_curve = fit(OF',R_e','cubicinterp');

end