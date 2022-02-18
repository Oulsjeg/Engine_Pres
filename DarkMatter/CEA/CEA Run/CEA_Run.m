% This routine runs the NASA CEA program for computing combustion products
% and properties for rocket engine applications. 
%
% DEPENDENCIES
%--------------------------------------------------------------
% This routine is executed from the Master.m main function file
%
% INPUTS
%--------------------------------------------------------------
% OF: Oxidizer to fuel ratio
% Pcc: Combustion chamber pressure
% SupAr: Supersonic area ratio
% PcPe: Chamber to exit pressure ratio
% fuel: Fuel data
% ox: Oxidizer data
%
% OUTPUTS
%--------------------------------------------------------------
% out: Output data
%--------------------------------------------------------------

function out = CEA_Run(OF, Pcc, SupAr, fuel, ox)

test = CEA;

test.OF       = OF;         %Set input O/F Ratio
test.pressure = Pcc;        %Set input Chamber Pressure
test.presUnit = 'psia';     %Set input pressure unit
test.supar    = SupAr;      % Supersonic Area Ratio

test.setFuel(fuel{1}, ...   % Set fuel
             fuel{2}, ...   % Set fuel percetnage
             fuel{3});      % Set fuel temperature
         
test.setOxid(ox{1}, ...     % Set oxidizer
             ox{2}, ...     % Set oxidizer percetnage
             ox{3});        % Set oxidizer temperature
        
ioinp = test.input.rocket;  % Set inputs
out = test.run;             % Run CEA

end
