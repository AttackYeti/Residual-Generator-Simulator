%% Settings
clear
close all

limits = 100:3500;   % In Seconds
deltat = 0.1;           % Size of timestep between readings (in seconds)
normalize = "true"; % Only enabled if [normalize = "true";]

% Define file names
filename = '2019-05-15 Fault Detection III DATABASE Supervisor.csv';
DAQfileName ='2019-05-15 Fault Detection III DATA.csv';
% Manually define specific fault initiation times
CTAH = [1638 1680 1690 1735 1760 1785];
HL = [1945 2036];

%% Start
disp("Beginning Simulation Script...");
tic

disp("Importing Data...");
DATA = importData(DAQfileName);

%% Heater parameter values
% Define specific heat capacities
disp("Defining Models...");
cpS = 500; %J/kgK 304 SS
cpF = 1800; %J/kgK 1518 + 2.82*T [deg C] Dowtherm A taking rough average
            % temp of 100 degrees Celsius
% Define fluid density
rho = 994.9; %kg/m3 (Dowtherm A) taking rough average temp of 100 degrees
             % Celsius and saturated liquid properties(from Dowtherm A Heat
             % Transfer Fluid Product Technical Data, Dow Chemical Company)

% Heater mass and geometry
mS = 5.9832; %kg, Outer shell mass multiplied by fraction in control
                  % volume
mT = (2.047 + 0.430)*0.741; %kg, Inner tube mass (new) + tape, multiplied
                            % by effective fraction in control volume
A = 2*pi*(1.5/2)*62.992*0.0254^2; % m2, Outer shell (inside) SA, 2 * pi * r * 
                              % L, then convert from inches to meters

% Fluid volume, mass and geometry
%mdot = 0.18; %kg/s, mass flow rate
V_S = pi*(1.5/2)^2*62.992*0.0254^3; % m3, Outer shell volume (with no inner tube)
V_T = ((255.8 + 56.1)/100^3)*0.741; %m3, Inner tube volume (new) tube + tape,
                                    % multipled by effective fraction in 
                                    % control volume
mF = (V_S - V_T)*rho; % kg
mI = mT + mF; % kg, mass of tube plus mass of fluid to get total inner mass

% Compute weighted average specific heat capacity
cpFav = (cpS*mT + cpF*mF)/mI; % J/kgK

% Heat transfer coefficient (guess)
h = 469; % Improved Heat Transfer and Volume Scaling through Novel Heater 
         % Design (Lukas, Kendrick, and Peterson) 
         
%% Heater model
modelDef.type = 'Symbolic';
modelDef.x = {'CX101','CX102','CX103','CX104','CX111',...
    'CX112','CX113','CX114','WT10','BT11','ST14E','ST14W','ST14N','PH',...
    'FM40','dTFH','tauH','TinH','TFH','TSH','dTSH','TFHbar','dTFHbar'};
modelDef.f = {'fH','fFM40'};
modelDef.z = {'uH','y1','y2','y3','y4','y5','y6','y7','y8','y9','y10',...
    'y11','y12','y13','y14'};
modelDef.parameters = {'hSH','ASH','MFH','cpFH','MSH','cpSH','VFH'};

% Define parameter values
params.hSH = h;
params.ASH = A;
params.MFH = mF;
params.cpFH = cpF;
params.MSH = mS;
params.cpSH = cpS;
params.VFH = V_S - V_T;
modelDef.parameter_values = params;

syms(modelDef.x{:})
syms(modelDef.f{:})
syms(modelDef.z{:})
syms(modelDef.parameters{:})

modelDef.rels = {y1==CX101,y2==CX102,y3==CX103,y4==CX104,... % e1,e2,e3,e4
    y5==CX111,y6==CX112,y7==CX113,y8==CX114,... %e5,e6,e7,e8
    y9==WT10,y10==BT11,y11==ST14E,y12==ST14W,... %e9,e10,e11,e12
    y13==ST14N, y14==FM40 + fFM40,uH==PH,... %e13,e14,e15
    dTFH==(hSH*ASH)/(MFH*cpFH)*(TSH-TFH),... %e16
    dTSH==(hSH*ASH)/(MSH*cpSH)*(TFH-TSH) + PH/(MSH*cpSH) + fH,... %e17
    tauH==VFH/FM40,... %e18
    DiffConstraint('dTFH','TFH'),... %e19
    DiffConstraint('dTSH','TSH'),... %e20
    TFH==(CX101+CX102+CX103+CX104+CX111+CX112+CX113+CX114)/8,... %e21
    TSH==(ST14E+ST14W+ST14N)/3,... %e22
    TinH==(WT10+BT11)/2,... %e23
    dTFHbar==(1/tauH)*(TinH - TFH),... %e24
    TFHbar==(TinH + TFH)/2,... %e25
    DiffConstraint('dTFHbar','TFHbar')... %e26
    };
HeaterModel = DiagnosisModel(modelDef);
HeaterModel.name = 'CIETheater';

% clear temporary variables from workspace
clear( modelDef.x{:} )
clear( modelDef.f{:} )
clear( modelDef.z{:} )
clear( modelDef.parameters{:} )

clear modelDef params

%% CTAH parameter values
% Define specific heat capacities
cpS = 390; %J/kgK Copper
cpF = 1771.8; %J/kgK 1518 + 2.82*T [deg C] Dowtherm A taking rough
                      % average temp of 90 degrees Celsius 
                      % (100 inlet, 80 outlet);
% Define fluid density
rho = 994.9; %kg/m3 (Dowtherm A) taking rough average temp of 100 degrees
             % Celsius and saturated liquid properties(from Dowtherm A Heat
             % Transfer Fluid Product Technical Data, Dow Chemical Company)
rho_S = 8960; %kg/m3 (Copper) 

% CTAH mass and geometry
iD = 0.466*0.0254; % inner diameter of piping, meters
oD = 0.5*0.0254; % outer diameter of piping, meters
iA = 0.17055*0.0254^2; % flow area, square meters
iL = 373.8*0.0254; % flow length, meters
iV = iA*iA; % area*length, cubic meters
iC = 2*pi*(iD/2); % inner circumference, 2*pi*d/2, meters
iSA = iC*iL; % inner surface area
sV = pi*((oD - iD)/2)^2*iL; % solid volume, pi*r^2*L
oC = 2*pi*(oD/2); % outer circumferece, 2*pi*od/s, meters
oSA = oC*iL; % outer surface area

% Heat transfer coefficient (guess)
hF = 292.15; % W/m2K, from Dane's notes and studies, calibrated at
             % 8kW power from 80-100 deg C

%% CTAH model
modelDef.type = 'Symbolic';
modelDef.x = {'WT40','BT41','WT42','BT43','AT01',...
    'dTFC','TinC','TFC','TSC','TFC','dTSC','TSC','TinftyC',...
    'FM40','hinftyC','tauC','omegaC','dTFCbar','TFCbar'};
modelDef.f = {'fC'};
modelDef.z = {'uC','y1','y2','y3','y4','y5'};
modelDef.parameters = {'hFC','AFC','MFC','cpFC','MSC','cpSC','AinftyC',...
    'VFC','omegaCfxn'};

% Define parameter values
params.hFC = hF;
params.AFC = iSA;
params.MFC = iV*rho;
params.cpFC = cpF;
params.MSC = sV*rho_S;
params.cpSC = cpS; 
params.AinftyC = oSA;
params.VFC = iV + 0.003; % Corrected for pipe elbows (transit time) 
params.omegaCfxn = 1122.3*(0.006295784*1.225); % from Jason's paper on
                                               % Model-based Non-Linear 
                                               % PID Control of a Heat 
                                               % Exchanger
modelDef.parameter_values = params;

syms(modelDef.x{:})
syms(modelDef.f{:})
syms(modelDef.z{:})
syms(modelDef.parameters{:})

modelDef.rels = {y1==WT40,y2==BT41,y3==WT42,y4==BT43,... % e1,e2,e3,e4
    y5==AT01,uC==omegaC,hinftyC==omegaCfxn*omegaC + fC,... %e5,e6,e7
    tauC==VFC/FM40,TFC==(BT41+WT40)/2,TinC==(BT43+WT42)/2,... %e8,e9,e10
    TinftyC==AT01,... %e11
    dTFC==(hFC*AFC)/(MFC*cpFC)*(TSC-TFC),... %e12
    dTSC==(hFC*AFC)/(MSC*cpSC)*(TFC-TSC)+(hinftyC*AinftyC)/(MSC*cpSC)*(TinftyC-TSC),... %e13
    dTFCbar==(1/tauC)*(TinC-TFC),... %e14
    TFCbar==(TinC+TFC)/2,... %e15
    DiffConstraint('dTFC','TFC'),... %e16
    DiffConstraint('dTSC','TSC'),... %e17
    DiffConstraint('dTFCbar','TFCbar')... %e18
    };
CTAHModel = DiagnosisModel(modelDef);
CTAHModel.name = 'CTAH';

% clear temporary variables from workspace
clear( modelDef.x{:} )
clear( modelDef.f{:} )
clear( modelDef.z{:} )
clear( modelDef.parameters{:} )

clear modelDef params

%% DRACS parameter values
cp = 1630.8; % [J/kgK] Specific heat capacity of Dowtherm A 1518 + 2.82*T [deg C] Dowtherm A taking rough average
           % temp of 40 degrees Celsius
rho = 1043.8; %kg/m3 (Dowtherm A) taking rough average temp of 40 degrees
             % Celsius and saturated liquid properties(from Dowtherm A Heat
             % Transfer Fluid Product Technical Data, Dow Chemical Company)
h = 3.66*16.2/(0.0139319*2); % [W/m2K] convection heat transfer coefficient DowA-to-pipe for circular tube, const. surf. temp, laminar, fully developed (4.36*k/D) for 304L SS 1" schedule 10
h_infty = 45; % [W/m2K] convection heat transfer coefficient air-to-pipe (from Dane's calculations in UCBTH-CIET-Nodalization - Dane de Wet Model - Final_20190311)
A1 = 2*pi*0.0139319; % [m2] inside pipe surf area (2*pi*r_i for 1" schedule 10 SS pipe)
A2 = 2*pi*0.0167005; % [m2] outside pipe surf area (2*pi*r_o for 1" schedule 10 SS pipe)
M = pi*0.0139319^2*1043.8; % [kg] mass per unit pipe length between WT and BT thermocouples of Dowtherm A - pi*r^2*rho
L = 3.800475; % [m] length of pipe between DHX outlet static mixer and TCHX inlet (RELAP nodalization (CIET_Model_Parameters_Structural_Mass_Dowtherm_Inertia))
Ltc = 0.0762; %[m] length between TC points

R = 1/(h*A1) + 1/(h_infty*A2); % [K/W] effective thermal resistance DowA to air (neglecting insulation)

%% DRACS model
modelDef.type = 'Symbolic';
modelDef.x = {'BT63','WT62','WT64','AT02',...
    'BT65','dT1','T1','dT2','T2','Tav',...
    'Tinfty'};
modelDef.f = {'fT1','fT2'};
modelDef.z = {'y1','y2','y3','y4','y5'};
modelDef.parameters = {'R','h','A','cp','M','Ltc'};

% Define parameter values
params.R = R;
params.M = M;
params.h = h;
params.A = A1*L;
params.cp = cp;
params.Ltc = Ltc;
modelDef.parameter_values = params;

syms(modelDef.x{:})
syms(modelDef.f{:})
syms(modelDef.z{:})
syms(modelDef.parameters{:})

modelDef.rels = {y1==BT63,y2==WT62,... % e1,e2
    y3==BT65,y4==WT64,y5==AT02,... %e3,e4,e5
    T1==(WT62+BT63)/2,... %e6
    Tinfty==AT02,... %e7
    T2==(WT64+BT65)/2,... %e8
    dT1==1/(R*cp*M*Ltc)*(Tinfty-T1)+fT1,... %e9
    dT2==1/(R*cp*M*Ltc)*(Tinfty-T2)+fT2,... %e10
    DiffConstraint('dT1','T1'),... %e11
    DiffConstraint('dT2','T2'),... %e12
    T2==T1+(h*A*R)*(Tinfty-Tav),... %e13
    };
DRACSModel = DiagnosisModel(modelDef);
DRACSModel.name = 'DRACS TC';

%% Hot Leg parameter values
cp = 1771.8; %J/kgK 1518 + 2.82*T [deg C] Dowtherm A taking rough
                      % average temp of 90 degrees Celsius 
                      % (100 inlet, 80 outlet);
rho = 1003.2; %kg/m3 (Dowtherm A) taking rough average temp of 90 degrees
             % Celsius and saturated liquid properties(from Dowtherm A Heat
             % Transfer Fluid Product Technical Data, Dow Chemical Company)
h = 3.66*16.2/(0.0139319*2); % [W/m2K] convection heat transfer coefficient DowA-to-pipe for circular tube, const. surf. temp, laminar, fully developed (4.36*k/D) for 304L SS 1" schedule 10
h_infty = 45; % [W/m2K] convection heat transfer coefficient air-to-pipe (from Dane's calculations in UCBTH-CIET-Nodalization - Dane de Wet Model - Final_20190311)
A1 = 2*pi*0.0139319; % [m2] inside pipe surf area (2*pi*r_i for 1" schedule 10 SS pipe)
A2 = 2*pi*0.0167005; % [m2] outside pipe surf area (2*pi*r_o for 1" schedule 10 SS pipe)
M = pi*0.0139319^2*1043.8; % [kg] mass for 3" pipe length between WT and BT thermocouples of Dowtherm A - pi*r^2*rho
L = 2.0217892; % [m] length of pipe between thermocouples (heater outlet to CTAH inlet, Dane's CIET nodalization (UCBTH-CIET-Nodalization - Dane de Wet Model - Final_20190311))
Ltc = 0.0762; %[m] length between TC points

R = 1/(h*A1) + 1/(h_infty*A2); % [K/W] effective thermal resistance DowA to air (neglecting insulation)

%% Hot Leg model
modelDef.type = 'Symbolic';
modelDef.x = {'BT12','WT13','BT43','WT42','AT02',...
    'dT1','T1','dT2','T2','Tav',...
    'Tinfty'};
modelDef.f = {'fT1','fT2'};
modelDef.z = {'y1','y2','y3','y4','y5'};
modelDef.parameters = {'R','h','A','cp','M','Ltc'};

% Define parameter values
params.R = R;
params.M = M;
params.h = h;
params.A = A1*L;
params.cp = cp;
params.Ltc = Ltc;
modelDef.parameter_values = params;

syms(modelDef.x{:})
syms(modelDef.f{:})
syms(modelDef.z{:})
syms(modelDef.parameters{:})

modelDef.rels = {y1==BT12,y2==WT13,... % e1,e2
    y3==BT43,y4==WT42,y5==AT02,... %e3,e4,e5
    T1==(WT13+BT12)/2,... %e6
    Tinfty==AT02,... %e7
    T2==(WT42+BT43)/2,... %e8
    dT1==1/(R*cp*M*Ltc)*(Tinfty-T1)+fT1,... %e9
    dT2==1/(R*cp*M*Ltc)*(Tinfty-T2)+fT2,... %e10
    DiffConstraint('dT1','T1'),... %e11
    DiffConstraint('dT2','T2'),... %e12
    T2==T1+(h*A*R)*(Tinfty-Tav),... %e13
    };
HotLegModel = DiagnosisModel(modelDef);
HotLegModel.name = 'HotLeg';
disp("Simulation Model Defined...");

%% Run Residual Generators
disp("Beginning Residual Generator Simulation...");

lastOutput = 0;
fprintf('Simulation Status:  ');
format short
statusInterval = 400;

for i = 1:height(DATA)
    timestep = DATA(i,:);
    
    % Code for Printing Progress to the Command Window
    if i - lastOutput > statusInterval
        lastOutput = i;
        percentage = (lastOutput/height(DATA))*100;
        if lastOutput < 2*statusInterval
            fprintf('%.2f %% ',percentage);
        elseif percentage < 10
            fprintf('\b\b\b\b\b\b\b%.2f %% ',percentage);
        else
            fprintf('\b\b\b\b\b\b\b\b%.2f %% ',percentage);
        end
    end
    
    if i == 1
        [CIETHeaterState.z,   heaterState_fb]                          = CIETHeaterStates(timestep, timestep);
        [CTAHState.z,         CTAHState_fb]                            = CTAHStates(timestep, timestep);
        [HotLegState.z,       heaterOutState_fb,   CTAHInState_fb]     = HotLegStates(timestep);
        [DRACSState.z,        DRACSOutState_fb,    TCHXInState_fb]     = DRACSStates(timestep); 
    else
        % Assign states
        CIETHeaterState.z                           = CIETHeaterStates(timestep, timestep);
        CTAHState.z                                 = CTAHStates(timestep, timestep);
        HotLegState.z                               = HotLegStates(timestep);
        DRACSState.z                                = DRACSStates(timestep);
        
        % Calculate residuals
        [res.rHeater(i),    heaterState_fb]         = ResGenCIETHeater_1_4( CIETHeaterState.z, heaterState_fb,       HeaterModel.parameter_values,   0.1);                
        [res.rCTAH(i),      CTAHState_fb]           = ResGenCTAH_1_5(       CTAHState.z,       CTAHState_fb,         CTAHModel.parameter_values,     0.1);
        [res.rHeaterOut(i), heaterOutState_fb]      = ResGenHeaterOut_2_1(  HotLegState.z,     heaterOutState_fb,    HotLegModel.parameter_values,   0.1);
        [res.rDRACSOut(i),  DRACSOutState_fb]       = ResGenDRACSOut_1_1(   DRACSState.z,      DRACSOutState_fb,     DRACSModel.parameter_values,    0.1);
        [res.rCTAHIn(i),    CTAHInState_fb]         = ResGenCTAHIn_1_1(     HotLegState.z,     CTAHInState_fb,       HotLegModel.parameter_values,   0.1);
        [res.rTCHXIn(i),    TCHXInState_fb]         = ResGenTCHXIn_2_1(     DRACSState.z,      TCHXInState_fb,       DRACSModel.parameter_values,    0.1);
        % Optional normalization
        if normalize == "true"
            % Normalize residuals
   
            res.rHeater(i)      = res.rHeater(i)      / CIETHeaterState.z(1);
            res.rCTAH(i)        = res.rCTAH(i)        / -1 * CTAHState.z(1) * CTAHModel.parameter_values.omegaCfxn;
            res.rHeaterOut(i)   = res.rHeaterOut(i)   / HotLegState.z(1);
            res.rDRACSOut(i)    = res.rDRACSOut(i)    / DRACSState.z(3);
            res.rCTAHIn(i)      = res.rCTAHIn(i)      / HotLegState.z(3);
            res.rTCHXIn(i)      = res.rTCHXIn(i)      / DRACSState.z(1);

        end     
    end
end

fprintf('\b\b\b\b\b\b\b\b%.1f %% \n',100);
disp("Fault simulation complete...");
timeElapsed = toc;
fprintf('Simulation completed in %.2f seconds...\n', timeElapsed);
fprintf("%d time steps analyzed...", height(DATA));

%% Begin Data Analysis Section
disp("Running Data Analysis...");

%% Import supervisor data

% Detect table import options
Topts = detectImportOptions(filename);

% Read data file into MATLAB table format
T = readtable(filename,Topts);

%% Define variable names

% Initialize string array of data column headings
DataName = strings([1,length(T.Properties.VariableNames)]);

% Extract string array of data column headings
for i = 1:length(T.Properties.VariableNames)
    DataName(i) = convertCharsToStrings(T.Properties.VariableNames{i});
end

%% Extract data from supervisor file

% Read DHX Temp Out Residual Generator
DHXTempOutRes = T{1:end,find(DataName == 'DHXTempOutResGen')};
Spoof = T{1:end,find(DataName == 'SpoofEnabled')};

%% Detect spoof changes
SpoofChanges = DATA.Time(find(abs((diff(Spoof))) == 1));
SpoofChanges(isnan(SpoofChanges)==1) = [];
SpoofChanges(3)=1253;

%% Redefine relevant variables for each residual generator

HeaterResSigs   = res.rHeater;
HeaterResTemps  = [DATA.CX101, DATA.CX102, DATA.CX103, DATA.CX104, ...
                   DATA.CX111, DATA.CX112, DATA.CX113, DATA.CX114, ...
                   DATA.WT10,  DATA.BT11,  DATA.ST14E, DATA.ST14W, ...
                   DATA.ST14N];
HeaterResOther  = [DATA.DesiredPower, DATA.PowerOut, DATA.FM40];

CTAHResSigs     = res.rCTAH; 
CTAHResTemps    = [DATA.WT40, DATA.BT41, DATA.WT42, DATA.BT43, DATA.AT01];
CTAHResOther    = [DATA.FM40, DATA.CTAHFreq];

HotLegResSigs   = [res.rHeaterOut; res.rCTAHIn]; 
HotLegResTemps  = [DATA.BT12, DATA.WT13,DATA.WT42, DATA.BT43, DATA.AT02];

DRACSResSigs    = [DHXTempOutRes, res.rTCHXIn]; 
DRACSResTemps   = [DATA.WT62, DATA.BT63, DATA.WT64, DATA.BT65, DATA.AT02];

Faults          = [CTAH, HL, SpoofChanges'];

%% Calculate Time Vector
len = height(DATA);
endTime = len*deltat;
timeScale = linspace(0,endTime,len);

%% Heater Faults
disp("Plotting Heater Faults...");
figure(1)
% Residual and faults
subplot(3,1,1)
plot(timeScale,HeaterResSigs)
hold on
ylabel('Residual Value')
yyaxis right
for i = 1:size(Faults,2)
    line([Faults(i) Faults(i)], ylim,'LineStyle','--');
    hold on
end
for i = 1:length(SpoofChanges)
    line([SpoofChanges(i) SpoofChanges(i)], ylim,'LineStyle','--');
    hold on
end
xlim([limits(1), limits(end)])
%ylim([min(HeaterResSigs(limits)),max(HeaterResSigs(limits))])
legend('Heater Fault Signal','Spoofs and Faults')
title('Heater Residual Generator')
set(gca,'YTick',[])
% Relevant temperatures
subplot(3,1,2)
for i = 1:size(HeaterResTemps,2)
    plot(timeScale, HeaterResTemps(:,i))
    hold on
end
xlim([limits(1), limits(end)])
title('Relevant Heater Temps')
ylabel('Temp [C]')
legend('CX101', 'CX102', 'CX103', 'CX104', 'CX111', 'CX112',...
    'CX113', 'CX114', 'WT10', 'BT11', 'ST14E', 'ST14W', 'ST14N')
% Relevant control signals
subplot(3,1,3)
yyaxis left
plot(timeScale, HeaterResOther(:,1))
ylabel('Power Out [W]')
hold on
plot(timeScale, HeaterResOther(:,2))
yyaxis right
for i = 3:size(HeaterResOther,2)
    plot(timeScale, HeaterResOther(:,i))
    hold on
end
xlim([limits(1), limits(end)])
ylabel('Flow Rate [kg/s]')
title('Relevant Heater Control Signals')
xlabel('Time [s]')
%ylim([min(res.rHeater) max(res.rHeater)])
legend('Desired Power','Output Power','Flow Rate')

%% CTAH Plots
disp("Plotting CTAH data...");
figure(2)
subplot(3,1,1)
plot(timeScale,CTAHResSigs);
hold on

title('CTAH Residual Generator')
ylabel('Residual Value')
xlim([limits(1), limits(end)])
%ylim([min(CTAHResSigs(limits)) max(CTAHResSigs(limits))])
for i = 1:size(Faults,2)
    line([Faults(i) Faults(i)], ylim,'LineStyle','--');
    hold on
end
for i = 1:length(SpoofChanges)
    line([SpoofChanges(i) SpoofChanges(i)], ylim,'LineStyle','--');
    hold on
end
legend('CTAH Fault Signal','Spoofs and Faults')
subplot(3,1,2)
for i = 1:size(CTAHResTemps,2)
    plot(timeScale, CTAHResTemps(:,i))
    hold on
end
xlim([limits(1), limits(end)])
title('Relevant CTAH Temps')
ylabel('Temp [C]')
legend('WT40', 'BT41', 'WT42', 'BT43', 'AT01', 'AT02')
subplot(3,1,3)
plot(timeScale, CTAHResOther(:,1))
hold on
ylabel('Flow Rate [kg/s]')
yyaxis right
plot(timeScale, CTAHResOther(:,2))
xlim([limits(1), limits(end)])
title('Relevant CTAH Control Signals')
ylabel('Frequency [Hz]')
xlabel('Time [s]')
legend('Primary Loop Flow Rate','CTAH Frequency')

%% Hot Leg Plots
disp("Plotting hot leg data...");
figure(3)
subplot(2,1,1)
plot(timeScale, HotLegResSigs(1,:))
plot(timeScale, HotLegResSigs(2,:))
hold on
for i = 1:size(Faults,2)
    line([Faults(i) Faults(i)], ylim,'LineStyle','--');
    hold on
end
for i = 1:length(SpoofChanges)
    line([SpoofChanges(i) SpoofChanges(i)], ylim,'LineStyle','--');
    hold on
end
xlim([limits(1), limits(end)])
%ylim([min(min((HotLegResSigs(limits)))), max(max(HotLegResSigs(limits)))])
legend('Heater Outlet Fault Signal','CTAH Inlet Fault Signal','Spoofs and Faults')
title('Hot Leg Residual Generator')
ylabel('Residual Value')
subplot(2,1,2)
for i = 1:size(HotLegResTemps,2)
    plot(timeScale, HotLegResTemps(:,i))
    hold on
end
xlim([limits(1), limits(end)])
title('Relevant Hot Leg Temps')
legend('BT12', 'WT13', 'WT42', 'BT43', 'AT01', 'AT02')
xlabel('Time [s]')
ylabel('Temp [C]')

%% DRACS Plots
disp("Plotting DRACS data...");
figure(4)
subplot(2,1,1)
plot(timeScale, DRACSResSigs)
hold on
for i = 1:size(Faults,2)
    line([Faults(i) Faults(i)], ylim,'LineStyle','--');
    hold on
end
for i = 1:length(SpoofChanges)
    line([SpoofChanges(i) SpoofChanges(i)], ylim,'LineStyle','--');
    hold on
end
xlim([limits(1), limits(end)])
%ylim([min(min((DRACSResSigs(limits)))), max(max(DRACSResSigs(limits)))])
legend('DHX Outlet Fault Signal','TCHX Inlet Fault Signal','Spoofs and Faults')
title('DRACS Residual Generator')
ylabel('Residual Value')
subplot(2,1,2)
for i = 1:size(DRACSResTemps,2)
    plot(timeScale, DRACSResTemps(:,i))
    hold on
end
xlim([limits(1), limits(end)])
legend('WT62', 'BT63', 'WT64', 'BT65', 'AT01', 'AT02')
title('Relevant DRACS Temps')
xlabel('Time [s]')
ylabel('Temp [C]')

%% End Script
finishTime = toc;
fprintf("Script completed in %.2f seconds",finishTime);









