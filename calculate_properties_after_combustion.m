function [T, gamma, cp] = calculate_properties_after_combustion(mDotN2O,mDotFuel,N2OT,fuelT)

%mDotN2O = 0.1163;
%mDotFuel = 0.0251;
%N2OT    = 288;
%fuelT = 2000;


%--------------------Chemical Data------------------------%
%Paraffin is C25H52 
 stoichiometricMolarOFRatio  = 76;
 molarWeightN2O              = 44;  %kg/kmol
 molarWeightParaffin         = 352; %kg/kmol
 molarWeightH2O              = 18;
 molarWeightCO2              = 44;
 molarWeightN2               = 14;
%The following values come from "CRC Handbook of Chemistry and Physics 77th
%Edition" (can be found in Magdalene library)

%Enthalpies of formation at standard temperature and pressure
 deltaHFN2O                 = 81.6;     %kJ/mol
 deltaHFH2O                 = -241.8;   %kJ/mol
 deltaHFCO2                 = -393.5;   %kJ/mol
 deltaHFN2                  = 0.0;      %not a typo
 
%Enthalpy of combustion
 deltaHCParaffin            = 16192;   %kJ/mol
 
%Enthalpy of fusion for paraffin
 deltaHFusParaffin          = 0.59566;  %J/mol
 
%Standard cp ( note - in the future better models could be used to take
%into account real gas properties)
 molarCpN2O                      = 38.6;     %J/molK
 molarCpH2O                      = 33.6;     %J/molK (assumes that all water in the chamber exists as vapour)
 molarCpCO2                      = 37.1;     %J/molK
 molarCpN2                       = 29.1;
 molarCpParaffin                 = 750;  %J/molK
 
%Standard molar cvs for the gases
 molarCvN2O                     = 20.4;     %J/molK
 molarCvCO2                     = 28.9;     %J/molK
 molarCvN2                      = 20.8;     %J/molK
 molarCvH2O                     = 26.28;    %J/molK

%We need to calculate lambda in order to check whether we have to do lean
%or rich calculations

 molarInletFlowN2O           = 1000 * mDotN2O/molarWeightN2O;       %mol/s
 molarInletFlowParaffin      = 1000 * mDotFuel/molarWeightParaffin; %mol/s
 
 %We declare the variables for the flow of each species that we expect
 %to see in the exit stream
 
 molarExitFlowN2O           = 0;
 molarExitFlowH2O           = 0;
 molarExitFlowCO2           = 0;
 molarExitFlowN2            = 0;
 molarExitFlowParaffin      = 0;
 
 %Standard Conditions
 standardTemp       = 298; %K
 
 lambda     = (molarInletFlowN2O/molarInletFlowParaffin)/(stoichiometricMolarOFRatio);
 
 if ( lambda >= 1 )
     %We are running lean
     molarExitFlowCO2       = molarInletFlowParaffin * 25;
     molarExitFlowH2O       = molarInletFlowParaffin * 26;
     molarExitFlowN2O       = molarInletFlowN2O - ( molarExitFlowH2O + 2* molarExitFlowCO2);
     molarExitFlowN2        = molarInletFlowN2O - molarExitFlowN2O;
 else
     %we are running rich
     molarExitFlowN2        = molarInletFlowN2O;
     molarExitFlowParaffin  = molarInletFlowParaffin - molarInletFlowN2O/76;
     molarExitFlowH2O       = molarInletFlowN2O + 50 * ( molarExitFlowParaffin - molarInletFlowParaffin);
     molarExitFlowCO2       = 25*(molarInletFlowParaffin-molarExitFlowParaffin);
 end
     
 
 %To simplify things we break it down into a heat transfer to get
 %everything to standard conditions ( and make the reasonable assumption of
 %negligible enthalpy variation with pressure ). We then react at standard
 %conditions. We then input the heat.
 
 %Initial heat transfer on nitrous
 qDotN2O       = molarInletFlowN2O * molarCpN2O * (N2OT - standardTemp);   %J/s
 
 qDotParaffin  = molarInletFlowParaffin * molarCpParaffin * (fuelT - standardTemp);  %J/s

 qDotStandardCond   = qDotN2O + qDotParaffin; %This is positive when cooling the inlet flows down
 
 %We perform the combustion
 
 paraffinBurned             = molarInletFlowParaffin - molarExitFlowParaffin;
 heatReleaseFromCombustion  = paraffinBurned * deltaHCParaffin*1000;
 
 %We calculate the energy put back into the exit gases
 qDotReheat     = qDotStandardCond + heatReleaseFromCombustion;
 
 %We calculate the properties of the exit gases (NEGLECTING THE PARAFFIN
 %CONTRIBUTION)
 
 totalMolarGas      = molarExitFlowCO2 + molarExitFlowH2O + molarExitFlowN2 + molarExitFlowN2O;
 molarExitFracN2    = molarExitFlowN2   /totalMolarGas;
 molarExitFracCO2   = molarExitFlowCO2  /totalMolarGas;
 molarExitFracN2O   = molarExitFlowN2O   /totalMolarGas;
 molarExitFracH2O   = molarExitFlowH2O  /totalMolarGas;
 
 molarCpExit        = molarExitFracN2*molarCpN2 + molarExitFracCO2*molarCpCO2 + molarExitFracH2O*molarCpH2O + molarExitFracN2O*molarCpN2O;
 molarCvExit        = molarExitFracN2*molarCvN2 + molarExitFracCO2*molarCvCO2 + molarExitFracH2O*molarCvH2O +  molarExitFracN2O*molarCvN2O; %J/molK
 molarWeightGas     = molarExitFracN2*molarWeightN2 + molarExitFracCO2*molarWeightCO2 + molarExitFracH2O*molarWeightH2O +  molarExitFracN2O*molarWeightN2O;
 gamma = molarCpExit/molarCvExit;
 cp = 1000*molarCpExit/(molarWeightGas);
 
 
 %We calculate the final temperature
 molarExitFlow  = (mDotN2O + mDotFuel)*1000/molarWeightGas;
 T = standardTemp + qDotReheat/(molarExitFlow*molarCpExit);

% gamma = 1.4;
% cp = 1005;
% T  = 2000;
 

end

