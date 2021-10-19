%% This code shows how the degassing model shown in Figure 7 was run for the dacitic composition.
% We use MELTS for Matlab version 1.0.
clear all 
close all
%% Load in path to MELTS for Matlab here (add your path in .... here).
addpath('c....\MELTS_Matlab-master\package')
warning('off', 'MATLAB:loadlibrary:cppoutput');

%%
% Setting final pressure in MPa
FP=0.15; % MPa
% final pressure in bars
FPbar=FP*10;


%% rhyolite-melts calculation for H2O
warning('off', 'MATLAB:loadlibrary:cppoutput');

% MELTS Dynamic 4 specifies rhyolite-MELTS 1.2.0, 
%It is the rhyolite-MELTS model with water properties completely replaced with the H2O-CO2 fluid 
%saturation model of Ghiorso and Gualda (2015)
liquidus = MELTSdynamic(4);

% This specifies the pressure to start the ascent at
pressure = 650;
% This specifies temperature, but it doesnt really matter, as you find 
% the liquidus anyway. Just needs to be above that. 
temperature = 1400.0;
%% Taking the Fe3/FeT ratio at this MgO from the MELTS model B
Fe3FeT=(0.3156*0.8998)/(0.3156*0.8998+7.7059)
FeT=(0.9159*0.8998+7.4187)
%% Setting bulk - average LL2 glass composition, Fe3FeT, CO2 and H2O from MELTS

bulk(1)=64.3225 %SiO2
bulk( 2) = 1.5491% TiO2
bulk( 3) =12.5764% Al2O3
bulk( 4) =FeT*Fe3FeT*1.1111% Fe2O3 
bulk( 5) = 0 % Cr2O3
bulk( 6) =FeT*(1-Fe3FeT)% FeO 
bulk( 7) =0.1610%  MnO
bulk( 8) = 1.0622%  MgO
bulk( 9) =0.0001; % NiO
bulk(10) =  0.0; %CoO
bulk(11) = 3.9603 %CaO
bulk(12) =3.7183 %Na2O
bulk(13) =2.4005 %K2O
bulk(14) = 0.5371 %P2O5
bulk(15)=1.78 % H2O wt%
bulk(16)=0.0241 % CO2 wt%
bulk(17) =  0.0;
bulk(18) =  0.0;
bulk(19) =  0.0;


%% This section findsthe liquidus
liquidus.engine.set('bulkComposition', bulk);
liquidus.engine.pressure = pressure;
liquidus.engine.temperature = temperature ;

% This means there is no heat loss from the system. 
liquidus.engine.setSystemProperties("Mode", "Isenthalpic");

% This bit finds the liquidus
liquidus.engine.findLiquidus;
temperature = liquidus.engine.temperature;

% makes a new melts dynamic. 
ptpath = liquidus.copyAndKeepOutput;
ptpath.engine.calcEquilibriumState;
disp(ptpath.engine.status.message);

% Adapted from Gleeson et al. in revision (GCA)
% While at more than 10* final presure, work down in 
% 10 bar steps, then as you aproach surface, work in 1 bar steps. 
while ptpath.engine.pressure >= FPbar
    if ptpath.engine.pressure>=10*FPbar
        ptpath = ptpath.addNodeAfter;
        ptpath.engine.pressure = ptpath.engine.pressure - 10;
        % First number determins model (2), 2 runs constant property. 
        % second number, 1 fractionates, 0 doesnt fractionate. 
        ptpath.engine.calcEquilibriumState(2, 1);
    end
    if ptpath.engine.pressure<10*FPbar
        ptpath = ptpath.addNodeAfter;
        ptpath.engine.pressure = ptpath.engine.pressure - 1;
        ptpath.engine.calcEquilibriumState(2, 0)
    end
end
%% Extracting variable useful things from the PT path


MgO= ptpath.getListProperty('dispComposition', 'liquid1', 'MgO');
H2O=ptpath.getListProperty('dispComposition', 'liquid1', 'H2O');
CO2=ptpath.getListProperty('dispComposition', 'liquid1', 'CO2');
FeO = ptpath.getListProperty('dispComposition', 'liquid1', 'FeO');
Fe2O3 = ptpath.getListProperty('dispComposition', 'liquid1', 'Fe2O3');
Al2O3 = ptpath.getListProperty('dispComposition', 'liquid1', 'Al2O3');
TiO2 = ptpath.getListProperty('dispComposition', 'liquid1', 'TiO2');
CaO = ptpath.getListProperty('dispComposition', 'liquid1', 'CaO');
SiO2= ptpath.getListProperty('dispComposition', 'liquid1', 'SiO2');
Na2O= ptpath.getListProperty('dispComposition', 'liquid1', 'Na2O');
K2O= ptpath.getListProperty('dispComposition', 'liquid1', 'K2O');
P2O5= ptpath.getListProperty('dispComposition', 'liquid1', 'P2O5');
MnO= ptpath.getListProperty('dispComposition', 'liquid1', 'MnO');
NiO=ptpath.getListProperty('dispComposition', 'liquid1', 'NiO');
FeOT=FeO+0.8998.*Fe2O3;
Fe3FeT=(2.*Fe2O3./159.69)./(FeO./71.844+2.*Fe2O3./159.69) 

P_1 = ptpath.getListProperty('pressure');
T_1 = ptpath.getListProperty('temperature');

%% This saves properties to feed into the Giordano Viscosity Calculator
% Here, we use Thermobar to perform calculations (Wieser et al. in prep,
% see github)
Giordano_Out=[SiO2' TiO2' Al2O3' FeOT' MnO' MgO' CaO' Na2O' K2O' P2O5' H2O' H2O'*0, T_1']
%% This saves properties provided in the supporting information. 
% some are duplicates of above. 
    
MgO_1 = ptpath.getListProperty('dispComposition', 'liquid1', 'MgO');
H2OM_1 = ptpath.getListProperty('dispComposition', 'liquid1', 'H2O');
CO2M_1 = ptpath.getListProperty('dispComposition', 'liquid1', 'CO2')';
% This gets the volume of the bulk system
Vol_bulk_1=ptpath.getListProperty('v','bulk')
mass_bulk_1 = ptpath.getListProperty('mass', 'bulk');
% This gives the mass and volume of the liquid
mass_melt_1 = ptpath.getListProperty('mass', 'liquid1');
Vol_melt_1=ptpath.getListProperty('v','liquid1')
% This gives the mass and volume of the fluid phase. 
Mass_vol_1=ptpath.getListProperty('mass','fluid1')
Vol_vol_1=ptpath.getListProperty('v','fluid1')'
% This replaces Nans with zeros.
Vol_vol_1(isnan( Vol_vol_1))=0
% This gets the MELTS viscosity (using the model of Shaw 1972
Vic=ptpath.getListProperty('viscosity')
% Convert viscosity to PaS
Vis_PaS=0.1*10.^(Vic)'

Cum_Vol_vol_1=ones(length(Vol_vol_1), 1)
Mass_plag_1=ptpath.getListProperty('mass','plagioclase1')
Mass_plag_1(isnan(Mass_plag_1))=0
%% Output for supplement
Excel_out=[P_1' CO2M_1 Vol_vol_1 Vol_melt_1' Vol_bulk_1' Vis_PaS]

%% H2O and CO2 vs P
figure
hold on
yyaxis left
plot(P_1, H2OM_1, '--g')
ylabel('H2O Melt, wt%')
yyaxis right
plot(P_1, CO2M_1, '-g')
xlabel('Pressure (bars)')
ylabel('CO2 Melt, wt%')
%% Volume % of fluid
Vol_per_Fluid=100.*(Vol_vol_1./Vol_bulk_1')
figure
hold on
plot(P_1, Vol_per_Fluid, '-g')
xlabel('Pressure (bars)')
ylabel('Volume % of fluid')
