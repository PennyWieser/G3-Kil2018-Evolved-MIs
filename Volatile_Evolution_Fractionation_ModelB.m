% This recreates Fig. 5 of Wieser et al. (2021). Note, that for some reason, with CO2=0.038, this
% crashes on my 16 GB ram PC, and has to be run on a 128 GB Ram (Thanks
% Matt Gleeson!)
% If it crashes, try changing the CO2 content so you can get the general
% gist of the code. 
clear all
close all
%% Load in your path to MELTS for MATLAB here (replace .. with path)
addpath('..\package')
warning('off', 'MATLAB:loadlibrary:cppoutput');

%% What type of melts, Pressure, temps, bulk H2O, XFe3
% MELTS Dynamic 4 specifies rhyolite-MELTS 1.2.0, 
%It is the rhyolite-MELTS model with water properties completely replaced with the H2O-CO2 fluid 
%saturation model of Ghiorso and Gualda (2015)
liquidus = MELTSdynamic(4);


display(liquidus.endMemberFormulas("bulk"))

pressure = 650; % in bars
temperature = 1400.0; % Celcius - but finds liquidus anyway
bulk(15) =  0.5; %H2O
bulk(16)=0.038 %CO2 content at saturation at these conditions

%% Fissure 8  starting composition

bulk(1)=50.6182% %SiO2
bulk( 2) = 2.5593%  % TiO2
bulk( 3) =13.3088%  % Al2O3
bulk( 4) =1.8223% % Fe2O3 
bulk( 5) = 0 % % Cr2O3
bulk( 6) =9.2939%; %FeO total
bulk( 7) =0.1828% MnO
bulk( 8) = 6.5894%  MgO
bulk( 9) =0; %%NiO 
bulk(10) =  0.0; %CoO
bulk(11) =10.7959%  %CaO
bulk(12) =2.3677% %%Na2O
bulk(13) =0.4759%  %K2O
bulk(14) = 0.2445%  %P2O5
bulk(17) =  0.0;
bulk(18) =  0.0;
bulk(19) =  0.0;

Tot=sum(bulk(1:19))

%% Find liquidus
liquidus.engine.set('bulkComposition', bulk); % Setting bulk composition as above
liquidus.engine.pressure = pressure; % Inputting pressure

liquidus.engine.temperature = temperature; % Inputting temperature

% Find liquidus
liquidus.engine.findLiquidus;
temperature = liquidus.engine.temperature;
disp(liquidus.engine.status.message);
display(temperature);


%% No fo2 buffer
liquidus.engine.setSystemProperties("Mode", "Fractionate Solids"); 
liquidus.engine.setSystemProperties("Mode", "Fractionate Fluids")

% Starts a new list (which reloads the library)
ptpath = liquidus.copyAndKeepOutput;

ptpath.engine.calcEquilibriumState; % How much of each mineral phase, what composition
disp(ptpath.engine.status.message);

% This sets up the mass of volatiles
massVolafter = 0;
massVolbefore = 0;

SumVol_mass=0
i=1;

ptpath2 = ptpath.copyAndKeepOutput

% We set up 2 paths, one before the system has cleared out (ptpath), one
% after  (ptpath2). Note, if ptpath is used, you end up with too much
% fractionated volatiles, as MELTS leaves some in the system to help with
% stability.
while ptpath.engine.temperature >= 726.85 +80+80 % End temperature in celcius
    i=i+1;
 
    ptpath = ptpath.addNodeAfter; % Putting point in for next temperature step
    ptpath.engine.temperature = ptpath.engine.temperature - 4; % Lowers tempearture by degrees
    ptpath2 = ptpath2.addNodeAfter; % Putting point in for next temperature step
    ptpath2.engine.temperature = ptpath2.engine.temperature - 4; % Lowers tempearture by degrees

    % Select 1 to get output after equilibration and before fractionation, 2 for output after fractionation
    % (either way, bulk composition will be updated after fractionation)
    ptpath.engine.calcEquilibriumState(1, 1);
    disp(ptpath.engine.status.message);
    massVolbefore(i) = ptpath.getNodeProperty([], 'mass','fluid1');
    volVolbefore(i) = ptpath.getNodeProperty([], 'v','fluid1');
      
   % This acounts for fluid left behind. E.g., output after fractionation. 
    ptpath2.engine.calcEquilibriumState(1, 2);
    massVolafter(i) = ptpath2.getNodeProperty([], 'mass','fluid1');
    volVolafter(i) = ptpath2.getNodeProperty([], 'v','fluid1');
    
    tf = isequal(ptpath.engine.status.message,'Quadratic iterations exceeded.');
    if tf ==1
        break
    end
    

        
end

%% Major element trajectories
MgO_melts4=ptpath.getListProperty('dispComposition', 'liquid1', 'MgO');
CaO_melts4=ptpath.getListProperty('dispComposition', 'liquid1', 'CaO');
FeO_melts4=ptpath.getListProperty('dispComposition', 'liquid1', 'FeO');
SiO2_melts4=ptpath.getListProperty('dispComposition', 'liquid1', 'SiO2');
TiO2_melts4=ptpath.getListProperty('dispComposition', 'liquid1', 'TiO2');
MnO_melts4=ptpath.getListProperty('dispComposition', 'liquid1', 'MnO');
Al2O3_melts4=ptpath.getListProperty('dispComposition', 'liquid1', 'Al2O3');
Na2O_melts4=ptpath.getListProperty('dispComposition', 'liquid1', 'Na2O');
K2O_melts4=ptpath.getListProperty('dispComposition', 'liquid1', 'K2O');
P2O5_melts4=ptpath.getListProperty('dispComposition', 'liquid1', 'P2O5');
Fe2O3_melts4=ptpath.getListProperty('dispComposition', 'liquid1', 'Fe2O3');
FeOt_melts4=FeO_melts4+0.89998*Fe2O3_melts4

%%  volume of the system at each step
vol_system=ptpath2.getListProperty('v','bulk')
%% Loss of volatiles
% This calculates the amount of volatiles actually fractionated in each step , e.g., volatiles
% before fractionation - volatiles after fractionation
Loss_Vol_mass=massVolbefore - massVolafter
Loss_Vol_vol=volVolbefore -volVolafter
% lets replace those pesky nans
Loss_Vol_vol(isnan(Loss_Vol_vol))=0
Loss_Vol_mass(isnan(Loss_Vol_mass))=0 
%% THis plot shows volatiles at start of step, and left after fractionation at each point. 
%You can see MELTS leaves a small mass in the system. we need to account
%for this, else you strip too much CO2 and H2O if you assume all the
%volatiles are fractionated.
figure
hold on
plot(MgO_melts4, massVolbefore, '.k')
plot(MgO_melts4, massVolafter, '.r')
xlabel('MgO Liq')
ylabel('Mass volatiles in system')
legend('Before frac', 'after frac')
%% Sum volatiles by mass

for i=1:length(Loss_Vol_mass)
    if i==1
    SumVol_mass(i)=Loss_Vol_mass(i)
    end
    if i>1
    SumVol_mass(i)=SumVol_mass(i-1)+Loss_Vol_mass(i)
    end
end
%% Lets check the lost volatiles sum is equal to the sum of the lost volatiles...
sumcheck=sum(Loss_Vol_mass)-SumVol_mass(end)
%% Same as above, but for volume not mass.
for i=1:length(Loss_Vol_mass)
    if i==1
    SumVol_vol(i)=Loss_Vol_vol(i)
    end
    if i>1
    SumVol_vol(i)=SumVol_vol(i-1)+Loss_Vol_vol(i)
    end
end
sumcheck=sum(Loss_Vol_vol)-SumVol_vol(end)
%% Volume of volatiles produced at each step (%) relative to volume of system
figure
plot(MgO_melts4, 100*Loss_Vol_vol./vol_system, '-k')
xlabel('MgO Melt')
ylabel('Vol volatiles lost / vol system at each step')
%% Total amount of volatiles produced up to each step relative to volume of system
figure
plot(MgO_melts4, 100*SumVol_vol./vol_system, '-k')
xlabel('MgO Melt')
ylabel('Tot Vol volatiles lost / vol remaining system (%)')
%% Total volume of volatiles produced relative to the initial volume of the system
figure
plot(MgO_melts4, 100*SumVol_vol./vol_system(1), '-k')
xlabel('MgO Melt')
ylabel('Tot Vol volatiles lost / initial vol system')

%% H2O sum in vapour
H2O_melts4=ptpath2.getListProperty('dispComposition', 'liquid1', 'H2O');
H2O_fluid=ptpath2.getListProperty('dispComposition', 'fluid1', 'H2O');
H2O_fluid(isnan(H2O_fluid))=0 %==0)

% Idea is that this foreloop calculates the total wt% lost to the fluid,
% e.g., for each step, wt% H2O lost = mass of fluid * wt% H2O in fluid/100.
% Add to sum of H2O lost  up to previous step.
for i=1:length(H2O_fluid)
    if i==1
    H2O_sum(i)=0
    end
    if i>1
    H2O_sum(i)=H2O_sum(i-1)+H2O_fluid(i)*Loss_Vol_mass(i)*0.01  % of H2O in fluid * Mass of fluid in that step so g of H2O
    end
end

%% Same Logic for CO2 sum in vapour
CO2_melts4=ptpath2.getListProperty('dispComposition', 'liquid1', 'CO2');
CO2_fluid=ptpath2.getListProperty('dispComposition', 'fluid1', 'CO2');
CO2_fluid(isnan(CO2_fluid))=0 %==0)


for i=1:length(CO2_fluid)
    if i==1
    CO2_sum(i)=0 
    end
    if i>1
    CO2_sum(i)=CO2_sum(i-1)+CO2_fluid(i)*Loss_Vol_mass(i)*0.01
    end
end
%% This is a good check, the sum of CO2 and H2O should be equal to the mass of volatiles loss

(H2O_sum(end)+CO2_sum(end))-sum(Loss_Vol_mass)

 %% Testing mass balance - Logic is that CO2* mass of volatiles + CO2* mass of liquid shoud equal CO2 put in initially
massMelt_melts4=ptpath2.getListProperty('mass','liquid1'); % Gets Mass of liquid phase at each step

CO2_initial=CO2_melts4(1) %0.0394
CO2_initial_g=CO2_melts4(1)*0.01*massMelt_melts4(1) % 0.0380
CO2_Final_Fluid=CO2_sum(end) %0.0566
CO2_Final_Melt_Conc=CO2_melts4(end) %0.0140
CO2_Final_Melt_totalg=CO2_melts4(end)*0.01*massMelt_melts4(end) %0.0027
Final_Sum=CO2_Final_Fluid+CO2_Final_Melt_totalg %0.0379 - E.g., good within rounding error

%% Same for H2O
H2O_initial=H2O_melts4(1) %0.5063
H2O_initial_g=H2O_melts4(1)*0.01*massMelt_melts4(1) % 0.5
H2O_Final_Fluid=H2O_sum(end) %0.0062
H2O_Final_Melt_Conc=H2O_melts4(end) %2.6023
H2O_Final_Melt_totalg=H2O_melts4(end)*0.01*massMelt_melts4(end) %0.4895
Final_SumH2O=H2O_Final_Fluid+H2O_Final_Melt_totalg %0.4957 - E.g., good within rounding error



%% Converting wt% in the fluid to mol%
CO2fluid_mol=CO2_fluid./44.01
H2Ofluid_mol=H2O_fluid./18.01
CO2fluid_molper=100.*CO2fluid_mol./(CO2fluid_mol+H2Ofluid_mol)
H2Ofluid_molper=100.*H2Ofluid_mol./(CO2fluid_mol+H2Ofluid_mol)
CO2fluid_wtperc=100.*CO2_fluid./(CO2_fluid+H2O_fluid)
H2Ofluid_wtperc=100.*H2O_fluid./(CO2_fluid+H2O_fluid)

%% Plot % by mass of CO2 and H2O lost to vpaour phase relative to initial mass of CO2 (e.g. conc in liquid * mass of liquid at start)
figure
hold on
box on
yyaxis left
plot(MgO_melts4, 100*CO2_sum/(CO2_melts4(1)*0.01*massMelt_melts4(1)), '-k')
ylabel('% CO2  loss')
yyaxis right
plot(MgO_melts4, 100*H2O_sum/(H2O_melts4(1)*0.01*massMelt_melts4(1)), ':k')
xlabel('MgO')
ylabel('% H2O  loss')
ax=gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%% Plot of mole proportion vs. MgO
figure
hold on
plot(MgO_melts4, CO2fluid_molper, '-r')
plot(MgO_melts4, CO2fluid_wtperc, ':r')

% yyaxis right
% plot(MgO_melts4, CO2fluid_molper, '-r')
% ylabel('mol % CO2 in fluid')
%% Total amount exsolved (mass). E.g., mass of CO2 over mass of melt in system at step 1
figure
hold on
yyaxis left
plot(MgO_melts4, 100*CO2_sum./massMelt_melts4(1), ':k')
plot(MgO_melts4, 100*H2O_sum./massMelt_melts4(1), '-k')
xlabel('MgO')
legend('CO2', 'H2O')
ylabel('100*mass volatile loss/initial mass of system') 
yyaxis right
plot(MgO_melts4, CO2fluid_molper, '-r')
ylabel('mol % CO2 in fluid')
%%
Pec_CO2_lost=100*CO2_sum/(CO2_melts4(1)*0.01*massMelt_melts4(1))
Pec_H2O_lost=100*H2O_sum/(H2O_melts4(1)*0.01*massMelt_melts4(1))
%% Viscosity using Shaw
viscosity=ptpath.getListProperty('viscosity');
Vis_PaS=0.1*10.^(viscosity)'
temp = ptpath.getListProperty('temperature');
temp_HT1987=20.1*MgO_melts4+1014
%% Output for Dataset_S1
% CaO in melt, MgO in Melt, Sum of CO2 exsolved so far, Sum of H2O exsolved
% so far, CO2 remainin in liquid, H2O remaining in liquid, mol perc of CO2
% in the vapour, mole perc of H2O in the vapour, mass of melt, 
Out=[temp' temp_HT1987' CaO_melts4' MgO_melts4' CO2_fluid' H2O_fluid' CO2_sum' H2O_sum' CO2_melts4' H2O_melts4' CO2fluid_molper' CO2fluid_wtperc' massMelt_melts4' Pec_CO2_lost' Pec_H2O_lost' Vis_PaS ]

%% Output for Thermobar Viscosity Code
Giordano_out=[SiO2_melts4' TiO2_melts4' Al2O3_melts4' FeOt_melts4' MnO_melts4' MgO_melts4' CaO_melts4' Na2O_melts4' K2O_melts4' P2O5_melts4' H2O_melts4']


