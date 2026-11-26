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
CaSOI = -3.2 %isnt this supposed to be 4.2 (according to line 3)


%% Load data (if txt file)
FullName        = fullfile('Data','ExampleDataSet.txt');
dataIn          = table2array(readtable(FullName));
[Nrows,Ncols]   = size(dataIn);                    % Determine size of array
NdatapointsperCycle = 720/0.2;                     % Nrows is a multitude of NdatapointsperCycle
Ncycles         = Nrows/NdatapointsperCycle;       % This must be an integer. If not checkwhat is going on
Ca              = reshape(dataIn(:,1),[],Ncycles); % Both p and Ca are now matrices of size (NCa,Ncycles)
p               = reshape(dataIn(:,2),[],Ncycles)*bara; 


%% Pressure preprocessing: pegging, averaging and filtering

ProcessedPressure = processMeasuredPressure(Ca, p, 1.00*bara);

% Individually pegged cycles for the original plots and work calculations
p = ProcessedPressure.PeggedPressure;

% Filtered mean pressure for pV and combustion analysis
pFiltered = ProcessedPressure.FilteredPressure;

%% Plotting 
f1=figure(1);
pp = plot(Ca,p/bara,'LineWidth',1);       
xlabel('Ca [deg]');ylabel('p [bar]');      
xlim([-360 360]);ylim([0.5 50]);                  
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
plot(V/dm^3,pFiltered/bara);
xlabel('V [dm^3]');ylabel('p [bar]');    
xlim([0 0.8]);ylim([0.5 50]);                     
set(gca,'XTick',[0:0.1:0.8],'XGrid','on','YGrid','on');  
title('Linear pV diagram')
addIdealDieselCycle(gca, Cyl);
totalWork = 0;

for cycle = 1:Ncycles
    W = trapz(V, p(:,cycle));
    totalWork = totalWork + W;
end

fprintf('Total work for %d cycles: %.2f J\n', Ncycles, totalWork);
fprintf('Average work per cycle: %.2f J\n', totalWork / Ncycles);


subplot(2,1,2)
loglog(V/dm^3,pFiltered/bara);
xlabel('V [dm^3]');ylabel('p [bar]'); 
xlim([0.02 0.8]);ylim([0 50]);                    
set(gca,'XTick',[0.02 0.05 0.1 0.2 0.5 0.8],...
    'YTick',[0.5 1 2 5 10 20 50],'XGrid','on','YGrid','on');   
title('Logarithmic pV diagram')

addIdealDieselCycle(gca, Cyl);

%% Combustion Analysis (aROHR, CA10, CA50, CA90)

T_IVC = 293.15; % [K] Intake temperature (we have to replace this with measured intake temperature)

% Initializing empty vectors
total_aROHR = zeros(NdatapointsperCycle, 1);
total_CA    = zeros(NdatapointsperCycle, 1);
total_gamma = zeros(NdatapointsperCycle, 1);

for cycle = 1:Ncycles
    p_cycle = smoothdata(p(:,cycle), 'sgolay', 21);
    CA_cycle = Ca(:,cycle);
    
    % Volume and derivatives for this cycle
    V_cycle = CylinderVolume(CA_cycle, Cyl);
    dVdCA   = gradient(V_cycle, CA_cycle);
    dpdCA   = gradient(p_cycle, CA_cycle);

    gamma = calcGamma(CA_cycle, p_cycle, V_cycle, CaIVC, CaSOI, T_IVC, SpS);
    
    % Apparent rate of heat release for this specific cycle
    aROHR_cycle = (gamma./(gamma-1)) .* p_cycle .* dVdCA + ...
                  (1./(gamma-1)) .* V_cycle .* dpdCA;
              
			  
    total_aROHR = total_aROHR + aROHR_cycle;
    total_CA    = total_CA + CA_cycle;
    total_gamma = total_gamma + gamma;
end

% Calculating average per cycle
aROHR  = total_aROHR / Ncycles;
CA_avg = total_CA / Ncycles;
gamma_avg = total_gamma / Ncycles;

% Restrict analysis to combustion region
combWindow = (CA_avg >= -20) & (CA_avg <= 80);

CA_comb    = CA_avg(combWindow);
aROHR_comb = aROHR(combWindow);

% Cumulative heat release
Qcum = cumtrapz(CA_comb,aROHR_comb);

% Normalize
Qnorm = Qcum ./ max(Qcum);

% Ensure monotonic interpolation
[Qunique,idx] = unique(Qnorm);

CA10 = interp1(Qunique,CA_comb(idx),0.10);
CA50 = interp1(Qunique,CA_comb(idx),0.50);
CA90 = interp1(Qunique,CA_comb(idx),0.90);

% Display results
fprintf('\nCombustion Analysis\n');
fprintf('CA10 = %.2f deg CA\n',CA10);
fprintf('CA50 = %.2f deg CA\n',CA50);
fprintf('CA90 = %.2f deg CA\n',CA90);

% Plot heat release
figure(3)

subplot(2,1,1)
plot(CA_cycle,aROHR,'LineWidth',1.5)
xlim([-45 135])  
ylim([-20 120])
xlabel('Crank Angle [deg]')
ylabel('aROHR [J/deg]')
title('Apparent Rate of Heat Release')
grid on

subplot(2,1,2)
plot(CA_comb,Qnorm,'LineWidth',1.5)
hold on
xline(CA10,'g','CA10')
xline(CA50,'b','CA50')
xline(CA90,'r','CA90')
xlabel('Crank Angle [deg]')
ylabel('Normalized Heat Release [-]')
title('Normalized Cumulative Heat Release')
legend('Heat Release','CA10','CA50','CA90')
grid on



%% Calculating KPI'S

% We need to change these to our own measurements
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

RunFuelAnalysis;

%% Combine original figures and fuel-analysis pages into one window

% Find the fuel-analysis window
mainWindow = findall(groot, ...
    'Type', 'figure', ...
    'Name', 'Theoretical fuel and injection-timing analysis');

assert(~isempty(mainWindow), ...
    'RunFuelAnalysis must run before this section.');

mainWindow = mainWindow(1);

% Find its existing tab group
tabGroups = findall(mainWindow, 'Type', 'uitabgroup');

assert(~isempty(tabGroups), ...
    'No tab group found in the fuel-analysis window.');

mainTabs = tabGroups(1);

pageNames = {
    'Measured cylinder pressure'
    'Measured pV diagrams'
    'Measured heat release'
    'Measured gamma'
};

for figureNumber = 1:4

    % Find an existing figure without creating a new one
    sourceFigure = findall(groot, ...
        'Type', 'figure', ...
        'Number', figureNumber);

    if isempty(sourceFigure)
        continue
    end

    sourceFigure = sourceFigure(1);

    % Never copy the destination window into itself
    if sourceFigure == mainWindow
        continue
    end

    % Copy axes, legends and colour bars together so their
    % relationships are preserved.
    children = get(sourceFigure, 'Children');

    keep = arrayfun(@(h) ...
        any(strcmp(get(h, 'Type'), ...
        {'axes', 'legend', 'colorbar'})), children);

    plotObjects = children(keep);

    if isempty(plotObjects)
        continue
    end

    newPage = uitab(mainTabs, ...
        'Title', pageNames{figureNumber});

    copyobj(plotObjects, newPage);

    % Close the old figure only after copying succeeds
    delete(sourceFigure);
end

set(mainWindow, ...
    'Name', 'Complete engine and fuel analysis', ...
    'NumberTitle', 'off');

drawnow;

%% Reorder pages

desiredOrder = {
    'Measured cylinder pressure'
    'Measured pV diagrams'
    'Measured heat release'
    'Fuel pV comparison'
    'Temperature and gamma'
    'SOI pV sweep'
    'Timing trends'
    'KPI results'
};

existingTabs = mainTabs.Children;
existingTitles = get(existingTabs, 'Title');

[found, indices] = ismember(desiredOrder, existingTitles);

assert(all(found), ...
    'Some page titles do not match the existing tabs.');

orderedTabs = existingTabs(indices);

% Keep any additional pages at the end
remainingTabs = existingTabs(~ismember(existingTabs, orderedTabs));

mainTabs.Children = [orderedTabs(:); remainingTabs(:)];

% Select the first page
mainTabs.SelectedTab = orderedTabs(1);

%% Add load-controlled analysis to the existing window
LoadStudyResults = faRunLoadStudy( ...
    Cyl, FuelAnalysisOptions, mainTabs);

LoadStudyKPI = LoadStudyResults.KPI;
disp(LoadStudyKPI);