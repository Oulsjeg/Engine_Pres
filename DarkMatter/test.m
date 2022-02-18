clear;
clc;
close all;
warning('off','all')                % Turn off warnings (not errors)

addfxnpath();

input           = loadInVar(1);

sim             = simulationClass(input);

sim.rocket.propulsion.propellants{1, 1}.propellantTank.getBlowdown();
sim.rocket.propulsion.propellants{2, 1}.propellantTank.getBlowdown();
