%% Info
% Toerental: 1500 RPM
% SOA van 4.2º voor TDC
% Resolutie van 0.2º CA
% Data voor 69 cycles (maximale van de Smetec, de OGO gensets kunnen in principe “onbeperkt” aan)
% 


%% init
clear all; clc;close all;
addpath( "Functions","Nasa");


%% Units
mm      = 1e-3;dm=0.1;
bara    = 1e5;
MJ      = 1e6;
kWhr    = 1000*3600;
volperc = 0.01; % Emissions are in volume percentages
ppm     = 1e-6; % Some are in ppm (also a volume- not a mass-fraction)
g       = 1e-3;
s       = 1;


%% Load NASA maybe you need it at some point?
% Global (for the Nasa database in case you wish to use it).
global Runiv
Runiv = 8.314;
[SpS,El]        = myload('Nasa\NasaThermalDatabase.mat',{'Diesel','O2','N2','CO2','H2O'});


%% Engine geom data
Cyl.Bore                = 104*mm;
Cyl.Stroke              = 85*mm;
Cyl.CompressionRatio    = 21.5;
Cyl.ConRod              = 136.5*mm;
Cyl.TDCangle            = 0;


% -- Valve closing events can sometimes be seen in fast oscillations in the pressure signal (due
% to the impact when the Valve hits its seat).
CaIVO = -355;
CaIVC = -135;
CaEVO = 149;
CaEVC = -344;
CaSOI = -3.2;


%% Load data (if txt file)
FullName        = fullfile('Data','ExampleDataSet.txt');
dataIn          = table2array(readtable(FullName));
[Nrows,Ncols]   = size(dataIn);                    % Determine size of array
NdatapointsperCycle = 720/0.2;                     % Nrows is a multitude of NdatapointsperCycle
Ncycles         = Nrows/NdatapointsperCycle;       % This must be an integer. If not checkwhat is going on
Ca              = reshape(dataIn(:,1),[],Ncycles); % Both p and Ca are now matrices of size (NCa,Ncycles)
p               = reshape(dataIn(:,2),[],Ncycles)*bara; 




%% Plotting 
f1=figure(1);
pp = plot(Ca,p/bara,'LineWidth',1);       
xlabel('Ca [deg]');ylabel('p [bar]');      
xlim([-360 360]);ylim([0 50]);                  
iselect = 10;                                    % Plot cycle 10 again in the same plot to emphasize it. Just to show how to access individual cycles.
line(Ca(:,iselect),p(:,iselect)/bara,'LineWidth',2,'Color','r');
YLIM = ylim;
% Add some extras to the plot
line([CaIVC CaIVC],YLIM,'LineWidth',1,'Color','b'); % Plot a vertical line at IVC. Just for reference not a particular reason.
line([CaEVO CaEVO],YLIM,'LineWidth',1,'Color','r'); % Plot a vertical line at EVO. Just for reference not a particular reason.
set(gca,'XTick',[-360:60:360],'XGrid','on','YGrid','on');  
title('All cycles in one plot.')


%% pV-diagram
V = CylinderVolume(Ca(:,iselect),Cyl);
f2 = figure(2);
subplot(2,1,1)
plot(V/dm^3,p(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]');    
xlim([0 0.8]);ylim([0.5 50]);                     
set(gca,'XTick',[0:0.1:0.8],'XGrid','on','YGrid','on');  
title('pV-diagram')

totalWork = 0;

for cycle = 1:Ncycles
    W = trapz(V, p(:,cycle));
    fprintf('Net work for cycle %d: %.2f J\n', cycle, W);
    totalWork = totalWork + W;
end

fprintf('Total work for %d cycles: %.2f J\n', Ncycles, totalWork);
fprintf('Average work per cycle: %.2f J\n', totalWork / Ncycles);


subplot(2,1,2)
loglog(V/dm^3,p(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]'); 
xlim([0.02 0.8]);ylim([0 50]);                    
set(gca,'XTick',[0.02 0.05 0.1 0.2 0.5 0.8],...
    'YTick',[0.5 1 2 5 10 20 50],'XGrid','on','YGrid','on');   
title('pV-diagram')

%% We need to change these to our own measurements
%-------------------------------------------------
RPM = 1500; % Rotations per minute [-]
Q_LHV = 43 * MJ; % Lower Heating Value of fuel [J/kg] 
mdot_fuel = 0.3e-3; % Fuel mass flow rate [kg/s]
mdot_CO2 = 0.9e-3; % CO2 mass flow rate [kg/s] 
mdot_NOx = 0.01e-3; % NOx mass flow rate [kg/s]
%-------------------------------------------------

%% KPI calculations
W_avg = totalWork / Ncycles; % Average net work per cycle [J]
f_cycle = (RPM / 60) / 2; % Cycles per second for a 4 stroke engine [Hz]
% Note that a 4 stroke engine takes 720 deg for a single cycle

P = W_avg * f_cycle; % power [W]

m_fuel = mdot_fuel / f_cycle; % Fuel mass per cycle [kg]
E_fuel = m_fuel * Q_LHV; % Fuel energy input per cycle [J]

Efficiency = W_avg / E_fuel; % Thermal efficiency [-]
bsfc = (mdot_fuel / P) / (g / kWhr); % Brake specific fuel consumption [g/kWhr]
bsCO2 = (mdot_CO2 / P) / (g / kWhr); % Brake specific CO2 emission [g/kWhr]
bsNOx = (mdot_NOx / P) / (g / kWhr); % Brake specific NOx emission [g/kWhr]

fprintf('Power: %.2f kW\n', P / 1000);
fprintf('Efficiency: %.1f %%\n', Efficiency * 100);
fprintf('bsfc: %.1f g/kWhr\n', bsfc);
fprintf('bsCO2: %.1f g/kWhr\n', bsCO2);
fprintf('bsNOx: %.1f g/kWhr\n', bsNOx);