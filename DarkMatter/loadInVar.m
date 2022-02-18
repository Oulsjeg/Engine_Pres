 function input   = loadInVar(rocket)
    switch rocket
        case 1
            % Engine Input Variables
            input.engine.name           = 'Bob';        % Engine name
            input.engine.Manu           = 'UTAT';       % Engine manufacturer
            input.engine.engfile        = 'enginename'; % Name of output engine file

            % Settings
            input.settings.efficiency   = 0.85;
            input.settings.flightType   = '2DOF';
            input.settings.cnv          = 6894.757;     % Psi to Pa
            input.settings.opt_alt_fac  = 0.55;         % Optimization altitude factor  [N/A]
            input.settings.alt_Bout     = 24300;        % Set burout altitude iterator      [m]
            input.settings.dPtol        = 200;
            input.settings.npr_inc      = 5e-7;
            input.settings.OF_i         = 1;            % initial OF value
            input.settings.OF_f         = 5;            % final OF value
            input.settings.num_OF       = 20;           % final OF timestep
            input.settings.OF           = 3;            % Design Ox/Fuel ratio      [double]
            input.settings.fluidModel   = 'coolprop';

            % Simulation parameters
            input.sim.numpt             = 200;          % Discretization points     [integer]
            input.sim.tBurn             = 7;            % Engine burn time          [s]
            input.sim.dt                = input.sim.tBurn/(input.sim.numpt + 1);
            input.sim.relax             = 0.3;          % Relaxation factor         [double, 0 < relax < ]
            input.sim.altConvCrit       = 50;           % Altitude convergence crit [m]
            input.sim.thetaL            = 2;            % Launch angle              [degrees]

            % Design Input Variables
            input.design.diameter       = 0.15;
            input.design.tankThickness  = 0.01;
            input.design.Pcc            = 350;          % Design chamber pressure   [psi]
            input.design.mDotox         = 1.4;          % Design Ox mdot            [kg/s]
            input.design.OF             = 3;
            input.design.Toxinit        = 278;          % Design Ox mdot            [kg/s]
            input.design.startPres      = 0;
            input.design.injCd          = 0.4;
            input.design.SupArInit      = 5;            % Supersonic Area Ratio initial guess

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % DEFINE PROPELLANT AND PRESSURANT TANKS
            % Fuel pressurant parameters
            input.fPres.name            = 'N2';         % Nitrogen                  [char]
            input.fPres.frac            = 100;          % Fraction                  [%]
            input.fPres.Tinit           = 278;          % Initial temperature       [K]
            input.fPres.rho             = 321;          % Density                   [kg/m^3]
            input.fPres.MW              = 28;           % Molar mass                [g/mol]
            input.fPres.Cp              = 0.28883e5;    % Heat capacity             [J/kmol K]
            input.fPres.mTank           = 2;            % READ FROM MASS BUDGET     [kg]
            input.fPres.lTank           = 1;            % READ FROM MASS MUDGET     [kg]
            input.fPres.vTank           = 0.02;         % READ FROM MASS BUDGET
            input.fPres.offset          = 0.2;          % Distance till next comp.  [m]
            input.fPres.mInit           = 3;
            input.fPres.Pinit           = 2500 * input.settings.cnv;

            % Fuel parameters
            input.fuel.isPropellant     = true;
            input.fuel.isPressurized    = true;
            input.fuel.pressurantOrder  = 'fwd';
            input.fuel.fluidtype        = 'Fuel';
            input.fuel.name             = 'C2H5OH';     % Ethanol
            input.fuel.pressurant       = 'fPres';
            input.fuel.blowdownMode     = 'constantMdot';
            input.fuel.frac             = 100;          % Fraction                  [%]
            input.fuel.Tinit            = 298;          % Initial temperature       [K]
            input.fuel.rho              = 789;          % Density                   [kg/m^3]
            input.fuel.MW               = 46.07;        % Molar mass                [g/mol]
            input.fuel.mTank            = 0;            % READ FROM MASS BUDGET
            input.fuel.lTank            = 0;            % READ FROM MASS MUDGET
            input.fuel.vTank            = 0.01;         % READ FROM MASS BUDGET
            input.fuel.order            = 2;            % Order inside rocket       [integer]
            input.fuel.offset           = 0.2;          % Distance till next comp.  [m]
            input.fuel.mInit            = 7;
            input.fuel.Pinit            = 525 * input.settings.cnv;


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Oxidizer pressurant parameters
            input.oxPres.name           = 'N2';         % Nitrogen
            input.oxPres.frac           = 100;          % Fraction                  [%]
            input.oxPres.Tinit          = 278;          % Initial temperature       [K]
            input.oxPres.rho            = 321;          % Density                   [kg/m^3]
            input.oxPres.MW             = 28;           % Molar mass                [g/mol]
            input.oxPres.Cp             = 0.28883e5;    % Heat capacity             [J/kmol K]
            input.oxPres.mTank          = 2;            % READ FROM MASS BUDGET
            input.oxPres.lTank          = 1;            % READ FROM MASS MUDGET
            input.oxPres.vTank          = 0.02;         % READ FROM MASS BUDGET
            input.oxPres.offset         = 0.2;          % Distance till next comp.  [m]
            input.oxPres.mInit          = 0.1;
            input.oxPres.Pinit          = 2500 * input.settings.cnv;

            % Oxidizer parameters
            input.ox.isPropellant       = true;
            input.ox.isPressurized      = true;
            input.ox.pressurantOrder    = 'fwd';
            input.ox.fluidtype          = 'Oxidizer';
            input.ox.name               = 'N2O';        % Nitrous Oxide
            input.ox.pressurant         = 'oxPres';
            input.ox.blowdownMode       = 'constantPressure';
            input.ox.frac               = 100;          % Fraction                  [%]
            input.ox.Tinit              = 278;          % Initial temperature       [K]
            input.ox.rho                = 881;          % Density                   [kg/m^3]
            input.ox.MW                 = 44.013;       % Molar mass                [g/mol]
            input.ox.mTank              = 10;           % READ FROM MASS BUDGET
            input.ox.lTank              = 2;            % READ FROM MASS BUDGET
            input.ox.vTank              = 0.02;         % READ FROM MASS BUDGET
            input.ox.order              = 1;            % Order inside rocket       [integer]
            input.ox.offset             = 0.2;          % Distance till next comp.  [m]
            input.ox.mInit              = 15.8;
            input.ox.Pinit              = 525 * input.settings.cnv;

    end

end
