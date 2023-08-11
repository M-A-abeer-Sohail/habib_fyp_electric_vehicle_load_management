%% Clear everything

clc;
close all;
clear;

%% Notes
% Save a dataset for each case

%   2 and 3 wheelers will have a different peak
%   They will have low charging powers and higher times
%   hence put up all night. They will also be more POPULAR.

% when dim = 1, remains the COLUMN, when dim = 2, remains the ROW

% VERY IMPORTANT: Cite the MATLAB codes
%   * Cite mult values from Satsangi

% IMPORTANT: for 20 pc no EV on bus 646

% power_13NodeTestFeeder -> Simulink command for IEEE 13 node

%% 20 pc set up

nptMult = 4;
npts = 24*nptMult;

titleFontSize = 20;
axisFontSize = 16;
lgndFontSize = 14;
plotWidth = 3;

% The details of making the new load file called Loadss.txt 

industkw = ceil([170 130 180 900 185*3 200 180 210]); % From the IEEE test case
industkvar = ceil(industkw.*[125 110 130 580 160*3 110 90 120]./industkw);

str_part=['New Load.645a Bus1=645.1.2   Phases=1 Conn=delta Model=2 kV=11 ',...
    'New Load.645b Bus1=645.2.3   Phases=1 Conn=delta Model=2 kV=11 ',...
    'New Load.645c Bus1=645.3.1   Phases=1 Conn=delta Model=2 kV=11 ',...
    "New Load.671 Bus1=671.1.2.3  Phases=3 Conn=Delta Model=1 kV=11 ",...
    'New Load.692  Bus1=692.1.2.3 Phases=3 Conn=Delta Model=5 kV=11  ',...
    'New Load.680a Bus1=680.1     Phases=1 Conn=Wye   Model=5 kV=6.35 ',...
    'New Load.680b Bus1=680.2     Phases=1 Conn=Wye   Model=5 kV=6.35 ',...
    'New Load.680c Bus1=680.3     Phases=1 Conn=Wye   Model=5 kV=6.35 ',...
    "New Load.634a Bus1=634.1     Phases=1 Conn=wye   Model=1 kV=0.23 ",...
    'New Load.634b Bus1=634.2     Phases=1 Conn=wye   Model=1 kV=0.23 ',...
    'New Load.634c Bus1=634.3     Phases=1 Conn=wye   Model=1 kV=0.23 ',...
    'New Load.646a Bus1=646.1     Phases=1 Conn=wye   Model=2 kV=0.23 ',...
    'New Load.646b Bus1=646.2     Phases=1 Conn=wye   Model=2 kV=0.23 ',...
    'New Load.652  Bus1=652.1     Phases=1 Conn=Wye   Model=1 kV=0.23 ',...
    'New Load.611  Bus1=611.3     Phases=1 Conn=Wye   Model=5 kV=0.23 ',...
    'New Load.675a Bus1=675.1     Phases=1 Conn=Wye   Model=1 kV=0.23 ',... 
    'New Load.675b Bus1=675.2     Phases=1 Conn=Wye   Model=1 kV=0.23 ',... 
    'New Load.675c Bus1=675.3     Phases=1 Conn=Wye   Model=1 kV=0.23 ']; 

housePowerFactor = 0.9; % Can change it later, we have assumed it

homes = ceil([60, 75, 55, 75, 95, 115, 117, 110, 120, 105]); 

% We are assuming the numbers of houses.
% Each house load is on a different single phase

% My goal is to have house load that has average load of 1.75 kW per
% household, AND peak load of 3.5 kW per household.
% Therefore, I took this value of a *rated house* load as
% 2.7 kW.

resLdRated = 2.7;

resLdNum = 10;
industLdNum = 8;

reskw = ceil(homes.*resLdRated); % Assumed from active power data from IoT dataset
reskvar = ceil(reskw*tan(acos(housePowerFactor)));

%%
close all;
noiseStd = 0.05; % For adding some randomness to test data

% the multiplier values for residential load, made using normpdf
% This peaks at 7pm and drops at close to 12 pm

multa = 1.15; % DIVIDES the houseMult values
multb = 1.25; % DIVIDES the industMult values

offseta = 0.2; % ADDS to houseMult values
offsetb = 0.3; % ADDS to industMult values

mults = load('multVals.mat');

% resMult = 5*normpdf(linspace(-0.1,0.5,npts),0.35,0.2) + noiseStd*5.*randn(1,npts);
% temp = load('resMult.mat'); % CITE FROM SATSANGI
% resMult = temp.ans; clearvars temp;
% resMult = resMult./max(resMult)/multa + offseta;
% resMult = resMult./max(resMult);

resMult = mults.resMult;
industMult = mults.industMult;

% the multiplier values for residential load, made using normpdf
% % This peaks at 10am and drops at close to 12 pm
% industMult = 5*normpdf(linspace(-0.7,0.5,npts),-0.1,0.3) + 5*noiseStd.*randn(1,npts);
% industMult = industMult./max(industMult)/multb + offsetb;
% industMult = industMult./max(industMult);

load busNames.mat;

basekVs = [132 11 11 0.4 11 0.4 11 0.4 0.4 11 0.4 11 11 11];

%Plotting the mult values
close all; plot(resMult); hold on; plot(industMult); hold off;
xticks(linspace(1,24*nptMult,7));
    xticklabels({'12:00 am','04:00 am','08:00 am','12:00 pm','04:00 pm','08:00 pm','12:00 am'});
    
reskwDaily = resMult'*reskw; % getting all 24*nptMult points for each houseLoad bus, kw
indkwDaily = industMult'*industkw; % getting all 24*nptMult points for each industLoad bus, kw

reskvarDaily = resMult'*reskvar; % getting all 24*nptMult points for each houseLoad bus, kvar
indkvarDaily = industMult'*industkvar; % getting all 24*nptMult points for each industLoad bus, kvar

resRatedkw = max(reskwDaily, [], 1);
industRatedkw = max(indkwDaily, [], 1);

resRatedkvar = max(reskvarDaily, [], 1);
industRatedkvar = max(indkvarDaily, [], 1);

% [houseRatedkw' houseRatedkvar'] 
% [industRatedkw' industRatedkvar']

makingLoads(str_part, [industkw reskw],...
    [industkvar reskvar]);

normalVals = load('normal.mat');

resLdkw = normalVals.resLdkw;
resLdkvar = normalVals.resLdkvar;
industLdkw = normalVals.industLdkw;
industLdkvar = normalVals.industLdkvar;

close all;
figure;
plot(resMult); hold on; plot(industMult); hold off;
figure;
plot(sum(reskwDaily, 2)); hold on;
plot(sum(indkwDaily, 2)); 
plot(sum(reskwDaily, 2) + sum(indkwDaily, 2)); hold off;

%% EV car numbers
% TABLE MAI ROUND OFF NHI DALI!!! -> Future Abeer: what?

pkCar = [224000 233000 244000];
khiCar = 0.2*pkCar;
evCarProp = [0.2 0.4 0.9];
evCarKhi = evCarProp.*khiCar;
evCar13 = [ceil(0.01*evCarKhi(1)) floor(0.01*evCarKhi(2:3))]; % This is the number of EV cars in the IEEE 13 bus
% evCar13 = [90 186 439];
 
pkBike = [2300000 2800000 3900000];
khiBike = 0.2*pkBike;
evBikeProp = [0.1 0.5 0.9];
evBikeKhi = evBikeProp.*khiBike;
evBike13 = floor(0.01*evBikeKhi); % This is the number of EV cars in the IEEE 13 bus
% evBike13 = [460 2800 7000];

%% EV load setup

numEVLds = 10;

evBusNames = [634, 675, 646, 652, 611]; % In an arbitrary but set order

Vi=ones(1,numEVLds); % Initial bus reference voltage, but I think I should
                     % lay a condition that relates to Haidar paper
                     % These are voltages per phase

% Will need to fix these...
evAtBus = readmatrix('data/evVeh20.csv');
carAtBus=evAtBus(:,1:end-2);
bikeAtBus=evAtBus(:,end-1:end);

evBikeV = [36 48]; % in volts
evBikeAh = [20 20]; % in Ah
evBikeChTime = [7 7]; % in hrs

% except EV volts, data taken from https://www.mobilityhouse.com/int_en/knowledge-center/charging-time-summary
% check battery voltages.xls
evCarV = [356 240 360 375 356 240 346 330 300 360 350]; % in Volts
evCarkwh = [64 28 27 30 64 4.4 8.8 16 12 24 30]; % in kwh
evCarChTime = [4.5 4.5 6 7.5 9.5 1.5 2.5 6 5 5.5 7]; % in hrs, wallbox

calcCarPavg(evCarkwh, evCarV, evCarChTime)
calcBikePavg(evBikeAh, evBikeV, evBikeChTime, 1)

evPowerFactor = 0.9;

busNumPh = [3 3 2 1 1];

str_part2 = ["New Load.EV_634a Bus1=634.1 Phases=1 Conn=Wye  Model=1 kV=0.23 vminpu=0.3 ",...
            'New Load.EV_634b Bus1=634.2 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
            'New Load.EV_634c Bus1=634.3 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
            'New Load.EV_675a Bus1=675.1 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
            'New Load.EV_675b Bus1=675.2 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
            'New Load.EV_675c Bus1=675.3 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
            "New Load.EV_646a Bus1=646.1 Phases=1 Conn=Wye  Model=1 kV=0.23 vminpu=0.3 ",...
            "New Load.EV_646b Bus1=646.2 Phases=1 Conn=Wye  Model=1 kV=0.23 vminpu=0.3 ",...
            "New Load.EV_652a Bus1=652.1 Phases=1 Conn=Wye  Model=1 kV=0.23 vminpu=0.3 ",...
            "New Load.EV_611c Bus1=611.3 Phases=1 Conn=Wye  Model=1 kV=0.23 vminpu=0.3 "];

% Generating EV loads for cars
[PevCar, QevCar] = getEvLoad(Vi, carAtBus, calcCarPavg(evCarkwh, evCarV, evCarChTime),...
    evPowerFactor, busNumPh, 1);
[PevBike, QevBike] = getEvLoad(Vi,bikeAtBus, calcBikePavg(evBikeAh, evBikeV, evBikeChTime, 1),...
    evPowerFactor, busNumPh, 1);

PevTot = PevCar + PevBike;
QevTot = QevCar + QevBike;

makingLoads([str_part str_part2], [industkw reskw PevTot'],...
        [industkvar reskvar QevTot']); % make the loadss.txt file

%%
close all;
% Setting up the load shape multipliers for EV charging
% carPdfMult = 1;
% 
% % I think the jaggedness in mult values resembles the
% % uncoordinated behavior of charging
% carpdf = normpdf(linspace(-0.8,0.5,npts),0,0.23) + 0.25*rand(1,npts); % this will multiply
%                                                                        % with the
%                                                                        % EV load to
%                                                                        % simulate gradual rise,
%                                                                        % also with some noise
% carpdf = +min(carpdf)*1.25 + carPdfMult*carpdf/max(carpdf);% normalizing, 
%                                             % and scaling 
%                                             % up/down with the mult
% carpdf = carPdfMult*carpdf/max(carpdf);
% 
% bikepdf = normpdf(linspace(-0.5,0.5,npts),0,0.27) + 0.1*rand(1,npts);% this will multiply
%                                                                        % with the
%                                                                        % EV load to
%                                                                        % simulate gradual rise,
%                                                                        % also with some noise
% bikepdf = -min(bikepdf)/2 + bikepdf/max(bikepdf);% normalizing, 
%                                % and scaling 
%                                % up/down with the mult
% bikepdf = bikepdf/max(bikepdf);
% 
% shiftA = 20;
% shiftB = 45;
% 
% carpdf = circshift(carpdf, shiftA);
% bikepdf = circshift(bikepdf, shiftB);

evMults = load('evMults.mat');
carpdf = evMults.carpdf;
bikepdf = evMults.bikepdf;

close all; plot(carpdf); hold on; % plot(circshift(carpdf, shiftA)); 
plot(bikepdf); % plot(circshift(bikepdf, shiftB)); 
hold off;
xticks(linspace(1,24*nptMult,7));
    xticklabels({'12:00 am','04:00 am','08:00 am','12:00 pm','04:00 pm','08:00 pm','12:00 am'});
    legend("1","2");
%%
% Used to set the voltage values in Vi for the new iteration when calculating EV load.
ViSetting = [[4,1];[4,2];[4,3];[6,1];[6,2];[6,3];[8,1];[9,3];[11,1];[11,2]]; 
                                                                       
% Declared once
lossesBus = {'Transformer.SUB';
            'Transformer.XFM1';
            'Line.650632';
            'Line.632671';
            'Line.671680';
            'Line.632633';
            'Line.632645';
            'Line.645646';
            'Line.692675';
            'Line.671684';
            'Line.684611';
            'Line.684652';
            'Line.671692'};

% Defined once
pNames = ["Transformer.SUB";'Transformer.XFM1';'Transformer.XFM2';'Transformer.XFM3';'Transformer.XFM4';
'Transformer.XFM5';'Line.650632';'Line.632671';'Line.671680';'Line.632633';'Line.632645';'Line.671684';
'Line.671692'; 'Capacitor.CAP1'; 'Capacitor.CAP2'; 'Load.645A'; 'Load.645B'; 'Load.645C'; 'Load.671';
'Load.692'; 'Load.680A'; 'Load.680B'; 'Load.680C'; 'Load.634A'; 'Load.634B'; 'Load.634C';
'Load.646A'; 'Load.646B'; 'Load.652'; 'Load.611'; 'Load.675A'; 'Load.675B';'Load.675C'];

ratedSubkw = 4788; % This is 20% of max load set as rated load.
ratedDisTr634kw = 617;
ratedDisTr652kw = 380;
ratedDisTr611kw = 374;
ratedDisTr675kw = 1086;
ratedDisTr646kw = 552;

%----------------------Junk Code---------------------

% DON'T DELETE - This shows the method of calculating those values
% Also, this was calculated from THE NORMAL CONDITIONS
% ratedLdkw = ceil(1.2*max(sum(reskwDaily, 2) + sum(indkwDaily, 2))); % 20% of max load we'll set to rated.
% ratedDisTr634kw = ceil(1.2*max(sum(reskwDaily(:,1:3), 2)));
% ratedDisTr652kw = ceil(1.2*max(reskwDaily(:,7)));
% ratedDisTr611kw = ceil(1.2*max(reskwDaily(:,6)));
% ratedDisTr675kw = ceil(1.2*max(sum(reskwDaily(:,8:end), 2)));
% ratedDisTr646kw = ceil(1.2*max(sum(reskwDaily(:,4:5), 2)));

%% 20pc condition all plots

% Arrays that'll hold per-phase voltage values
% These will be cleared for each case
V1pu = [];
V2pu = [];
V3pu = [];

% power flows for substation and distribution transformer
% These will be cleared for each case
subkw = [];
subkvar = [];

xfmkw652 = [];
xfmkvar652 = [];

xfmkw634 = [];
xfmkvar634 = [];

xfmkw611 = [];
xfmkvar611 = [];

xfmkw675 = [];
xfmkvar675 = [];

xfmkw646 = [];
xfmkvar646 = [];

subI = [];
xfm634I = [];
xfm675I = [];
xfm652I = [];
xfm611I = [];
xfm646I = [];
bus680I = [];
bus692I = [];

PevTot = [];
QevTot = [];

PLoss = [];
QLoss = [];

% This will be defined once as well
lossesAll = [];

LoadStart = 16;
% Used to set the loads phase wise in final matrix
phaseWiseVec = [4 5 6 0 0 1 2 3 1 2 3 1 2 1 3 1 2 3 1 2 3 1 2 3 1 2 1 3];
                                                                          
numTrafo = 6; % 5 distr. trafo's and one substation trafo

% clearvars reskwDaily reskvarDaily indkwDaily indkvarDaily

for i = 1:npts
%   residential house load
%     kwLoopRes = housekw*houseMult(i); kvarLoopRes = housekvar*houseMult(i); 
%     The particular power values for the i-th hour
%   industrial load
%     kwLoopInd = industkw*industMult(i); kvarLoopInd = industkvar*industMult(i);

    %residential house load
    kwLoopRes = reskwDaily(i,:); kvarLoopRes = reskvarDaily(i,:); % The particular power
    %industrial load                                              % values for the i-th hour
    kwLoopInd = indkwDaily(i,:); kvarLoopInd = indkvarDaily(i,:);
    
    % Generating EV loads
    [PevCarTemp, QevCarTemp] = getEvLoad(Vi, carAtBus, calcCarPavg(evCarkwh, evCarV, evCarChTime),...
        evPowerFactor, busNumPh, carpdf(i));
    
    [PevBikeTemp, QevBikeTemp] = getEvLoad(Vi, bikeAtBus, calcBikePavg(evBikeAh, evBikeV, evBikeChTime, 1),...
        evPowerFactor, busNumPh, bikepdf(i));
    
    PevTotTemp = PevCarTemp + PevBikeTemp;
    QevTotTemp = QevCarTemp + QevBikeTemp;
    
    makingLoads([str_part str_part2], [kwLoopInd kwLoopRes PevTotTemp'],...
        [kvarLoopInd kvarLoopRes QevTotTemp']); % make the loadss.txt file
    DailyLoadFlow; % compute the daily load flow from the files in OpenDSS
    V1pu = [V1pu; Vpu1hr(:,1)']; % Get the per-phase voltage values for phase a
    V2pu = [V2pu; Vpu1hr(:,2)']; % Get the per-phase voltage values for phase b
    V3pu = [V3pu; Vpu1hr(:,3)']; % Get the per-phase voltage values for phase c
    
    % Read the powers PER PHASE from the file
    [tempTrafokw, tempTrafokvar, tempLoadkw, tempLoadkvar]...
         = getPs(xlsread('MasterIEEE13_EXP_P_BYPHASE.CSV',1,'E2:J44'),...
         LoadStart, phaseWiseVec, numTrafo); % Extract the 
                                             % power values for all loads
                                             % and trafos
                                             
                                             
    PevTot = [PevTot; sum(tempLoadkw(end-numEVLds:end,:), 1)];
    QevTot = [QevTot; sum(tempLoadkvar(end-numEVLds:end,:), 1)];                                         

%     Again, method of calculation shown here
%     industLdkw = [industLdkw; sum(tempLoadkw(1:industLdNum,:), 1)];
%     industLdkvar = [industLdkvar; sum(tempLoadkvar(1:industLdNum,:), 1)];
%     
%     resLdkw = [resLdkw; sum(tempLoadkw(industLdNum+1:end,:), 1)];
%     resLdkvar = [resLdkvar; sum(tempLoadkvar(industLdNum+1:end,:), 1)];
    
    tempI = getCurr3ph1hr(xlsread('MasterIEEE13_EXP_CURRENTS.CSV',1,'L2:Q8'), 2, 7);
    tempBus680I = getCurr1ph1hr(xlsread('MasterIEEE13_EXP_CURRENTS.CSV',1,'B23:E25'));
    tempBus692I = getCurr3ph1hr(xlsread('MasterIEEE13_EXP_CURRENTS.CSV',1,'L2:Q15'), 14,14);
    
    subI = [subI; tempI(1,:)];
    xfm634I = [xfm634I; tempI(2,:)];
    xfm675I = [xfm675I; tempI(3,:)];
    xfm652I = [xfm652I; tempI(4,:)];
    xfm611I = [xfm611I; tempI(5,:)];
    xfm646I = [xfm646I; tempI(6,:)];
    bus680I = [bus680I; tempBus680I];
    bus692I = [bus692I; tempBus692I]; % Need this
    
    % Substation and distr. trafo.
    subkw = [subkw; tempTrafokw(1,:)]; % Store the kW power flow for substation
    subkvar = [subkvar; tempTrafokvar(1,:)]; % Store the kVAR power flow for substation
    
    xfmkw634 = [xfmkw634; tempTrafokw(2,:)]; % Store the kW power flow for dist. trafo. of 634
    xfmkvar634 = [xfmkvar634; tempTrafokvar(2,:)]; % Store the kVAR power flow for dist. trafo. of 634
    
    xfmkw675 = [xfmkw675; tempTrafokw(3,:)]; % Store the kW power flow for dist. trafo. of 675
    xfmkvar675 = [xfmkvar675; tempTrafokvar(3,:)]; % Store the kVAR power flow for dist. trafo. of 675
    
    xfmkw652 = [xfmkw652; tempTrafokw(4,:)]; % Store the kW power flow for dist. trafo. of 652
    xfmkvar652 = [xfmkvar652; tempTrafokvar(4,:)]; % Store the kVAR power flow for dist. trafo. of 652
        
    xfmkw611 = [xfmkw611; tempTrafokw(5,:)]; % Store the kW power flow for dist. trafo. of 611
    xfmkvar611 = [xfmkvar611; tempTrafokvar(5,:)]; % Store the kVAR power flow for dist. trafo. of 611
    
    xfmkw646 = [xfmkw646; tempTrafokw(6,:)]; % Store the kW power flow for dist. trafo. of 646
    xfmkvar646 = [xfmkvar646; tempTrafokvar(6,:)]; % Store the kVAR power flow for dist. trafo. of 646
    
    Vi = getNewVi(Vpu1hr, ViSetting, numEVLds);
    
%   Getting all loses for the i-th HOUR
                                                                                     
    [PLossTemp, QLossTemp] = getLinePower(readtable("MasterIEEE13_EXP_LOSSES.CSV")); 
    PLoss = [PLoss PLossTemp];
    QLoss = [QLoss QLossTemp];
    
    % clearvars -except subkw subkvar xfmkw xfmkvar str_part lds_kw lds_kvar normal_mult overld20_mult overld40_mult V1pu V2pu V3pu
end

save('ev20Pent','V1pu','V2pu','V3pu',...
    'xfmkw634', 'xfmkvar634', 'xfmkw675', 'xfmkvar675', 'xfmkw652', 'xfmkvar652',...
    'xfmkw611', 'xfmkvar611', 'xfmkw646', 'xfmkvar646',...
    'subkw','subkvar', 'resLdkw','industLdkw',...
    'resLdkvar','industLdkvar','lossesBus','PevTot','QevTot',...
    'resRatedkw','resRatedkvar','industRatedkw','industRatedkvar',...
    'subI','xfm634I','xfm646I','xfm675I','xfm652I','xfm611I','bus680I','bus692I',...
    'PLoss','QLoss');

clc;

%%
close all;

% Just testing here
% figure;
% plot(sum(subkw, 2),'Color','red');
% load('normal.mat');

V1pu(V1pu==0) = NaN;
V2pu(V2pu==0) = NaN;
V3pu(V3pu==0) = NaN;

%%
% Total load (residential + industrial) graphs
figure('Renderer', 'painters', 'Position', [200 -100 1100 700]);
hold on; 
plot(sum(resLdkw, 2) + sum(industLdkw, 2),'Color','blue','Marker','o','MarkerFaceColor','white',...
    'linewidth',plotWidth);
plot(sum(PevTot, 2) + sum(resLdkw, 2) + sum(industLdkw, 2),'Color','Red','Marker','o','MarkerFaceColor','white',...
    'linewidth',plotWidth);
hold off;

xticks(linspace(1,24*nptMult,7));
xticklabels({'12:00 am','04:00 am','08:00 am','12:00 pm','04:00 pm','08:00 pm','12:00 am'});
title('Daily load profile for all loads at 20% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
xlim([1 24*nptMult]); grid on;
ylim([(floor(min(sum(resLdkw, 2) + sum(industLdkw, 2))/500)*450),...
    (ceil(max(sum(PevTot, 2) + sum(resLdkw, 2) + sum(industLdkw, 2))/500)*500)])
lgnd = legend(["Non EV load","Total load plus EVs"],'Location','southoutside','orientation','horizontal');
    lgnd.FontSize = lgndFontSize; clear str;

busVSet = [2 4 5 6 8 14];  % Only the bus voltages I need

%% Voltages (check busNames for voltage values)
for j = 1:length(busVSet) % ignoring SOURCEBUS
    i = busVSet(j);
    figure('Renderer', 'painters', 'Position', [200 -100 1100 700])
    hold on;
    if ~isnan(V1pu(1,i))
         plot(V1pu(:,i),'linewidth',plotWidth);
    end
    if ~isnan(V2pu(1,i))
         plot(V2pu(:,i),'linewidth',plotWidth);
    end
    if ~isnan(V3pu(1,i))
         plot(V3pu(:,i),'linewidth',plotWidth);
    end
    plot(1.05*ones(1,npts),'linewidth',plotWidth,'linestyle','--','color','black');
    plot(0.95*ones(1,npts),'linewidth',plotWidth,'linestyle','--','color','black');
    hold off
    
    legndArr = {};
    if ~isnan(V1pu(1,i))
        legndArr(end+1) = {"Phase a"};
    end
    if ~isnan(V2pu(1,i))
         legndArr(end+1) = {"Phase b"};
    end
    if ~isnan(V3pu(1,i))
         legndArr(end+1) = {"Phase c"};
    end
    str1 = "1.05 p.u." + newline + "limit";
    str2 = "0.95 p.u." + newline + "limit";
    legndArr(end+1) = {str1};
    legndArr(end+1) = {str2};
    clearvars str1 str2
    lgnd = legend(string(legndArr),'Location','southoutside','orientation','horizontal');
    lgnd.FontSize = lgndFontSize; clear str;
    grid on;
    
    title(strcat("Voltage variation for bus", " ", busNames(i),", ","at 20 percent EV penetration (base kV_{\phi} ="," ",...
        num2str(round(basekVs(i)/sqrt(3),3))," ","kV_{rms})")); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Voltage in PU', 'FontSize', axisFontSize);
    ax = gca; ax.FontSize = axisFontSize;
    
    xticks(linspace(1,24*nptMult,7));
    xticklabels({'12:00 am','04:00 am','08:00 am','12:00 pm','04:00 pm','08:00 pm','12:00 am'});
    ylim([min([V1pu(:,i)' V2pu(:,i)' V3pu(:,i)' 0.95])-0.01 1.1])
    xlim([1 nptMult*24]);
end

%%
% Plotting 24-hr flow at substation
figure('Renderer', 'painters', 'Position', [200 -100 1100 700])
plot(subkw(:,1),'linewidth',plotWidth); hold on;
plot(subkw(:,2),'linewidth',plotWidth+0.2); 
plot(subkw(:,3),'linewidth',plotWidth+0.4);
plot(sum(subkw, 2),'linewidth',plotWidth+0.6);
plot(ones(1,npts)*ratedSubkw,'linewidth',plotWidth);
hold off;

xticks(linspace(1,24*nptMult,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Substation (bus 650) active power flow at 20% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Sum of" + newline + "all phases";
lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = lgndFontSize; clearvars str;
ylim([floor(min(subkw, [], 'all'))/500 ceil(max(sum(subkw, 2))/500)*500*1.25])
xlim([1 24*nptMult]);

%%
% Plotting 24-hr flow at dist. trafo. 634
figure('Renderer', 'painters', 'Position', [200 -100 1100 700])
plot(xfmkw634(:,1),'linewidth',plotWidth); hold on;
plot(xfmkw634(:,2),'linewidth',plotWidth+0.2); 
plot(xfmkw634(:,3),'linewidth',plotWidth+0.4); 
plot(sum(xfmkw634, 2),'linewidth',plotWidth+0.6);
plot(ones(1,npts)*ratedDisTr634kw,'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*nptMult,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Distribution transformer (634) power flow at 20% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Sum of" + newline + "all phases";
ylim([floor(min(min(xfmkw634), [], 'all')/500) ceil(max(sum(xfmkw634, 2)))*1.25])
xlim([1 24*nptMult]);
lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = lgndFontSize; clearvars str;

%%
% Plotting 24-hr flow at dist. trafo. 652
figure('Renderer', 'painters', 'Position', [200 -100 1100 700])
plot(xfmkw652(:,1),'linewidth',plotWidth*3); hold on;
plot(xfmkw652(:,2),'linewidth',plotWidth*1.5); 
plot(xfmkw652(:,3),'linewidth',plotWidth*0.8); 
plot(sum(xfmkw652, 2),'linewidth',plotWidth*1.15);
plot(ones(1,npts)*ratedDisTr652kw,'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*nptMult,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Distribution transformer (652) power flow at 20% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Sum of" + newline + "all phases";
ylim([floor(min((xfmkw652), [], 'all'))/500 ceil(max(sum(xfmkw652, 2)))*1.25])
xlim([1 24*nptMult]);
lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = lgndFontSize; clearvars str;

%%
% Plotting 24-hr flow at dist. trafo. 611
figure('Renderer', 'painters', 'Position', [200 -100 1100 700])
plot(xfmkw611(:,1),'linewidth',plotWidth); hold on;
plot(xfmkw611(:,2),'linewidth',plotWidth*1.12); 
plot(xfmkw611(:,3),'linewidth',plotWidth*1.14); 
plot(sum(xfmkw611, 2),'linewidth',plotWidth*1.16);
plot(ones(1,npts)*ratedDisTr611kw,'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*nptMult,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Distribution transformer (611) power flow at 20% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Sum of" + newline + "all phases";
ylim([floor(min((xfmkw611), [], 'all'))/500 ceil(max(sum(xfmkw611, 2)))*1.25])
xlim([1 24*nptMult]);
lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = lgndFontSize; clearvars str;

%%
% Plotting 24-hr flow at dist. trafo. 675
figure('Renderer', 'painters', 'Position', [200 100 1100 700])
plot(xfmkw675(:,1),'linewidth',plotWidth); hold on;
plot(xfmkw675(:,2),'linewidth',plotWidth+0.2); 
plot(xfmkw675(:,3),'linewidth',plotWidth+0.4); 
plot(sum(xfmkw675, 2),'linewidth',plotWidth+0.6);
plot(ones(1,npts)*ratedDisTr675kw,'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*nptMult,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Distribution transformer (675) power flow at 20% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Sum of" + newline + "all phases";
ylim([floor(min((xfmkw675), [], 'all'))/500 ceil(max(sum(xfmkw675, 2)))*1.25])
xlim([1 24*nptMult]);
lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = lgndFontSize; clearvars str;

%%
% Plotting 24-hr flow at dist. trafo. 646
figure('Renderer', 'painters', 'Position', [200 100 1100 700])
plot(xfmkw646(:,1),'linewidth',plotWidth); hold on;
plot(xfmkw646(:,2),'linewidth',plotWidth*1.12); 
plot(xfmkw646(:,3),'linewidth',plotWidth*1.14); 
plot(sum(xfmkw646, 2),'linewidth',plotWidth*1.16);
plot(ones(1,npts)*ratedDisTr646kw,'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*nptMult,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Distribution transformer (646) power flow at 20% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Sum of" + newline + "all phases";
ylim([floor(min((xfmkw646), [], 'all'))/500 ceil(max([sum(xfmkw646, 2); ratedDisTr646kw]))*1.15])
xlim([1 24*nptMult]);
lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = lgndFontSize; clearvars str;

%%
% Plotting 24-hr current flow at 2ndary terminal for sub
figure('Renderer', 'painters', 'Position', [200 100 1100 700])
plot(subI(:,1),'linewidth',plotWidth); hold on;
plot(subI(:,2),'linewidth',plotWidth*1.12); 
plot(subI(:,3),'linewidth',plotWidth*1.14); 
% plot(sum(xfmkw646, 2),'linewidth',plotWidth*1.16);
plot(ones(1,npts)*(1.1*5500/sqrt(3)/11),'linewidth',plotWidth);

% Rated value for "Dove" conductor should not be used with a trafo -> Future Abeer: why?

% plot(ones(1,npts)*(5500*1.5/11/sqrt(3)),'linewidth',plotWidth); -> What was this for?
hold off;
xticks(linspace(1,24*nptMult,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Substation (Sub) current flow at 20% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Current (A)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Normal Ampere rating";

% Max rating are from max current, normamps and emergamps
ylim([floor(min((subI), [], 'all'))/500 ceil(max([max(subI,[] ,'all'); 1.1*5500/sqrt(3)/11]))*1.15]) 
xlim([1 24*nptMult]);
lgnd = legend({'Phase a','Phase b','Phase c',"Normal Ampere rating"},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = lgndFontSize; clearvars str;

%%
% Plotting 24-hr current flow at 2ndary terminal for 634
figure('Renderer', 'painters', 'Position', [200 100 1100 700])
plot(xfm634I(:,1),'linewidth',plotWidth); hold on;
plot(xfm634I(:,2),'linewidth',plotWidth*1.12); 
plot(xfm634I(:,3),'linewidth',plotWidth*1.14); 
% plot(sum(xfmkw646, 2),'linewidth',plotWidth*1.16);
plot(ones(1,npts)*1000*1.1/0.4/sqrt(3),'linewidth',plotWidth);
% plot(ones(1,npts)*(1000*1.5/0.4/sqrt(3)),'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*nptMult,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Distribution transformer (634) current flow at 20% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Current (A)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Normal Ampere rating";

% Max rating are from max current, normamps
ylim([floor(min((xfm634I), [], 'all'))/500 ceil(max([max(xfm634I,[] ,'all'); 1000*1.1/0.4/sqrt(3)]))*1.15])  
xlim([1 24*nptMult]);
lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = lgndFontSize; clearvars str;

%%
% Plotting 24-hr current flow at 2ndary terminal for 675
figure('Renderer', 'painters', 'Position', [200 100 1100 700])
plot(xfm675I(:,1),'linewidth',plotWidth); hold on;
plot(xfm675I(:,2),'linewidth',plotWidth*1.12); 
plot(xfm675I(:,3),'linewidth',plotWidth*1.14); 
% plot(sum(xfmkw646, 2),'linewidth',plotWidth*1.16);
plot(ones(1,npts)*1000*1.1/0.4/sqrt(3),'linewidth',plotWidth);
% plot(ones(1,npts)*(1000*1.5/0.4/sqrt(3)),'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*nptMult,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Distribution transformer (675) current flow at 20% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Current (A)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Normal Ampere rating";

% Max rating are from max current, normamps
ylim([floor(min((xfm675I), [], 'all'))/500 ceil(max([max(xfm675I,[] ,'all'); 1000*1.1/0.4/sqrt(3)]))*1.15])
xlim([1 24*nptMult]);
lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = lgndFontSize; clearvars str;

%%
% Plotting 24-hr current flow at 2ndary terminal for bus 680
% IMPORTANT! I'll mention the current rating separately
figure('Renderer', 'painters', 'Position', [200 100 1100 700])
plot(bus680I(:,1),'linewidth',plotWidth); hold on;
plot(bus680I(:,2),'linewidth',plotWidth*1.12); 
plot(bus680I(:,3),'linewidth',plotWidth*1.14); 
% plot(ones(1,npts)*726,'linewidth',plotWidth); % CURRENT RATING!!! -> Why'd I remove this?
% plot(ones(1,npts)*(1000*1.5/0.4/sqrt(3)),'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*nptMult,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Bus 680 current flow at 20% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Current (A)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Normal Ampere rating";

% Max rating are from max current, normamps
ylim([floor(min((bus680I), [], 'all'))/500 ceil(max([max(bus680I,[] ,'all'); 0]))*1.15])  
xlim([1 24*nptMult]);
lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = lgndFontSize; clearvars str;

%%
% Plotting 24-hr current flow at 2ndary terminal for bus 692
% IMPORTANT! I'll mention the current rating separately
figure('Renderer', 'painters', 'Position', [200 100 1100 700])
plot(bus692I(:,1),'linewidth',plotWidth); hold on;
plot(bus692I(:,2),'linewidth',plotWidth*1.12); 
plot(bus692I(:,3),'linewidth',plotWidth*1.14); 
% plot(ones(1,npts)*726,'linewidth',plotWidth); % CURRENT RATING!!!
% plot(ones(1,npts)*(1000*1.5/0.4/sqrt(3)),'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*nptMult,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Bus 692 current flow at 20% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Current (A)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Normal Ampere rating";
ylim([floor(min((bus692I), [], 'all'))/500 ceil(max([max(bus692I,[] ,'all'); 0]))*1.15]) % Max rating are from max current,
                                                                                                      % normamps 
xlim([1 24*nptMult]);
lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = lgndFontSize; clearvars str;

%%
% Plotting 24-hr current flow at 2ndary terminal for bus 652
figure('Renderer', 'painters', 'Position', [200 100 1100 700])
plot(xfm652I(:,1),'linewidth',plotWidth); hold on;
plot(xfm652I(:,2),'linewidth',plotWidth*1.12); 
plot(xfm652I(:,3),'linewidth',plotWidth*1.14); 
plot(ones(1,npts)*1000*1.1/0.4/sqrt(3),'linewidth',plotWidth);
% We'll just mention it!!!
% plot(ones(1,npts)*(1000*1.5/0.4/sqrt(3)),'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*nptMult,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Distribution transformer (652) current flow at 20% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Current (A)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Normal Ampere rating";
ylim([floor(min((xfm652I), [], 'all'))/500 ceil(max([max(xfm652I,[] ,'all'); 0]))*1.15]) % Max rating are from max current,
                                                                                                      % normamps 
xlim([1 24*nptMult]);
lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = lgndFontSize; clearvars str;


%%
% This will contain all the losses
lossesAll = [lossesAll (sum(subkw, 2) - sum(resLdkw, 2) - sum(industLdkw, 2))];

% lossesPercent = [lossesPercent sum(lossesAll(:,1))*100/sum(sum(subkw, 2))];