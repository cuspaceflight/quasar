%Quasar.m
%Models Quasar using a quasi-equilibrium 0D model

%Written by Jack Brewster (jb803@cam.ac.uk) 2014

%Every time step we iteratively solve for chamber conditions assuming that
%the nozzle is sufficiently small to allow choking

%-------------------Nitrous Supply Data------------------------------%
supplyN2OPressure   = 50E5;     %50 bar
supplyN2ODensity    = 835.3;    %kg/m^3
supplyN2OTemp       = 288;      %K

%-------------------Injector Data------------------------------------%
injectorDiameter    = 1.9394E-3;   %m
injectorCD          = 0.6;      %dimensionless

%-------------------Combustion Chamber Data--------------------------%
initialPortDiameter = 25E-3;    %m
length              = 0.5;      %m
chamberPressure     = 35E5; %25 bar ( this is an initial pressure) - causing it to start at 50 bar will be bad as it cocks up  the combustion code
chamberTemperature  = 1500;     %K ( this is an initial temperature)
chamberDiameter     = 90E-3;    %m (sets a limit to port diameter growth) 


%-------------------Fuel Data----------------------------------------%
%Note that changing this will require changing of the following functions:
fuelDensity         = 900;      %kg/m^3
fuelRegressionRate  = 2E-3;     %m/s ( note that this may be overriden in later functions)
fuelOD              = 50E-3;    %m
fuelLength          = 0.3;      %m
fuelRefRegression   = 3E-3;     %m/s reference regression rate for plotting

%-------------------Nozzle Data--------------------------------------%
nozThroatDiameter   = 7.3E-3;     %m
nozExitDiameter     = 10E-3;

%-------------------Chemical Data------------------------------------%

%-------------------Run Data-----------------------------------------%
simEndTime          = 1e-3;        %Time to run till (s)
simTimeStep         = 1e-3;     %Time step ( time step may be dynamically altered in later versions) s
simCurrT            = 0;
simConvergenceLimit = 0.000001;     %When the change in the convergence variable ( pressure in this case ) is less than simConvergenceLim * convergenceVariable then we have converged
simMaxIterations    = 1000000;      %Maximum number of iterations before we give up and claim it's diverged
simDivergeFlag      = 0;        %A flag to indicate whether we've diverged
simMaxConvergeJump  = 0.7;      %The maximum fraction of the convergence variable we're allowed to change

%--------------------Variables for Results---------------------------%
resultPressure      = [];
resultTemperature   = [];
resultPortDiameter  = [];
resultTime          = [];

%------------------Initial Calculations------------------------------%
injectorArea        = 0.25*pi*injectorDiameter^2;
nozThroatArea       = 0.25*pi*nozThroatDiameter^2;
nozExitArea         = 0.25*pi*nozExitDiameter^2;
portDiameter        = initialPortDiameter;

%******************Main Loop*****************************************%
while simCurrT < simEndTime

    deltaPressureFraction       = 1;       %Forces at least a single loop of iteration
    iterationCounts             = 0;
    
    convergenceHistory  = chamberPressure;
    regressionRateHistory = fuelRegressionRate;
    
    while deltaPressureFraction > simConvergenceLimit
        
    
        %We firstly get the oxidizer mass flow
        mDotN2O     = injectorCD * injectorArea * sqrt( 2 * (supplyN2OPressure - chamberPressure ) * supplyN2ODensity);
    
        %We next get fuel mass flow
        %TODO: change this to use actual correlations
        %We may potentially need to iterate on this as well ( if there is a
        %temperature link)
        %Regression rate formula provided by slide 17 of https://web.stanford.edu/dept/aeroastro/cgi-bin/events/50thAnniversary/media/Karabeyoglu.pdf
        massFlux    = 1000*mDotN2O/(0.25*pi*(portDiameter*100)^2);
        fuelRegressionRate  = (0.488 * (massFlux^0.28))/1000;
        mDotF       = fuelRegressionRate * fuelDensity * fuelLength * pi * portDiameter;
    
        %We now perform the combustion analysis.
        %This is done via a seperate function in order to keep this code clean
        %and ensure that there is only one area which needs changing when
        %modifying the fuel
    
        [chamberTemperature,chamberGamma,chamberCp]     = calculate_properties_after_combustion(mDotN2O,mDotF,supplyN2OTemp,298);
    
        %We calculate the dimensionless mass flow coefficient wrt to choking
       
        
        dimensionlessMassFlow   = (chamberGamma/(sqrt(chamberGamma-1)))*((1+0.5*(chamberGamma-1))^(-0.5*((chamberGamma+1)/(chamberGamma-1))));
    
        mDotReq     = dimensionlessMassFlow*nozThroatArea*chamberPressure/(sqrt(chamberCp*chamberTemperature));
    
        mDotN2OReq  = mDotReq-mDotF;
    
        if mDotN2OReq<0
            disp('Backflow detected');
            simDivergeFlag = 1;
            break
        end
    

        chamberPressureReq      = supplyN2OPressure -0.5*(1/supplyN2ODensity)*(mDotN2OReq/(injectorCD*injectorArea))^2;
        
        %To prevent massive oscillatory jumps we force the maximum delta to be less than a fraction of the chamber pressure 
        deltaChamberPressure    = chamberPressureReq - chamberPressure;
        deltaPressureFraction   = abs(deltaChamberPressure/chamberPressure);

        
        %We update the chamber pressure for the new loop
        chamberPressure = chamberPressure+deltaChamberPressure; 
        
        convergenceHistory = [convergenceHistory chamberPressureReq];
        regressionRateHistory   = [regressionRateHistory fuelRegressionRate];
        

        
        iterationCounts     = iterationCounts + 1;
        if iterationCounts >= simMaxIterations
            simDivergeFlag = 1;
            disp('Chamber pressure failed to converge');
            
            break
        end
            
        
        
    end
    plot(1:(iterationCounts+1),convergenceHistory./supplyN2OPressure,'r-',1:(iterationCounts+1),regressionRateHistory./fuelRefRegression,'b-');
    %disp('iter done')
    
    %We check whether we diverged or suffered blow back
    if simDivergeFlag == 1
        disp('Diverged - ending sim')
        break
    end
    
    %We expand the port diameter
    portDiameter = portDiameter + simTimeStep * fuelRegressionRate * 2;
    
    if ( portDiameter > chamberDiameter )
        disp('Burnout');
        break
    end
    
    %We record all the data for plotting after simulation
    resultPressure      = [resultPressure chamberPressure];
    resultTemperature   = [resultTemperature chamberTemperature];
    resultPortDiameter  = [resultPortDiameter portDiameter];
    resultTime          = [resultTime simCurrT];
    
    simCurrT    = simCurrT + simTimeStep;
    
    
end


%plot(resultTime,resultPressure);


