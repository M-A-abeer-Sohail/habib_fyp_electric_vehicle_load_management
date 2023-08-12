%% Notes
% Save a dataset for each case 2 and 3 wheelers will have a different peak
% They will have low charging powers and higher times hence put up all night. They will also be more POPULAR.
% When dim = 1, remains the COLUMN, when dim = 2, remains the ROW
% VERY IMPORTANT: Cite the MATLAB codes
% power_13NodeTestFeeder -> Simulink command for IEEE 13 node
% When you use i as a count variable, use j for the complex number

%% Refactor notes
% Have meaningful variable names
% Remove duplication
% Get in working condition
% Get FYP results
% (Optional) Make it general purpose (this was our ambition)

%% Design notes
% Have only two files: normal condition and with EV
% with EV will take arguments for EV data

%% Normal condition set up

clc;
close all;
clear;

pointsMultiplier = 4;
nPoints = 24*pointsMultiplier;

titleFontSize = 20;
axisFontSize = 16;
legendFontSize = 14;
plotWidth = 3;

% The details of making the new load file called Loadss.txt 

industrialKwFromTestCase = ceil([170 130 180 900 185*3 200 180 210]); % From the IEEE test case
% TODO: Beautify this confusing calculation
industrialKvarFromTestCase = ceil(industrialKwFromTestCase.*[125 110 130 580 160*3 110 90 120]./industrialKwFromTestCase);

% Model set to PQ for all, because we are setting the kW

loadFileComponents=['New Load.645a Bus1=645.1.2   Phases=1 Conn=delta Model=2 kV=11 ',...
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

homePowerFactor = 0.9; 

nHomes = ceil([60, 75, 55, 75, 95, 115, 117, 110, 120, 105]);
% We are assuming the numbers of houses, and we actually worked backwards
% from the load we want and how much is the rated residential load.

% My goal is to have house load that has average load of 1.75 kW per
% household, AND peak load of 3.5 kW per household.
% Therefore, I took this value of a *rated house* load as
% 2.7 kW.

residentialRatedLoadKw = 2.7;

residentialKwLoad = ceil(nHomes.*residentialRatedLoadKw); % Assumed from active power data from IoT dataset
residentialKvarLoad = ceil(residentialKwLoad*tan(acos(homePowerFactor)));

nResidentialLoads = length(residentialKwLoad);
nindustrialLoads = length(industrialKwFromTestCase);

% fcn to make the Loadss.txt file for all hours
makingLoads(loadFileComponents, [industrialKwFromTestCase residentialKwLoad],...
    [industrialKvarFromTestCase residentialKvarLoad]);

if exist('multiplierValues.mat','var')
    multiplierValues = load('multiplierValues.mat');
    residentialMuls = multiplierValues.residentialMuls;
    industrialMults = multiplierValues.industrialMults;
else
    whiteNoise = 0.05; % For adding some randomness to test data

    % the multiplier values for residential load, made using normpdf
    % This peaks at 7pm and drops at close to 12 pm    
    multA = 1.15; % DIVIDES the houseMult values
    multB = 1.25; % DIVIDES the industMult values
    
    offsetA = 0.2; % ADDS to houseMult values
    offsetB = 0.3; % ADDS to industMult values

    % Generating the multipliers: residential loadshape multipliers
    temp = load('residentialMultsSatsangi.mat'); % CITED FROM SATSANGI
    residentialMuls = temp.ans; clearvars temp;
    residentialMuls = residentialMuls./max(residentialMuls)/multA + offsetA;
    residentialMuls = residentialMuls./max(residentialMuls);
    
    % Generating the multipliers: industrial loadshape multipliers
    % This peaks at 10am and drops at close to 12 pm
    industrialMults = 5*normpdf(linspace(-0.7,0.5,nPoints),-0.1,0.3) + 5*whiteNoise.*randn(1,nPoints);
    industrialMults = industrialMults./max(industrialMults)/multB + offsetB;
    industrialMults = industrialMults./max(industrialMults);
    
    % Saving those multipliers for later use
    save('multiplierValues','residentialMuls','industrialMults');
end

load busNames.mat;
perUnitBasekVs = [132 11 11 0.4 11 0.4 11 0.4 0.4 11 0.4 11 11 11];

%%
close all; plot(residentialMuls); hold on; plot(industrialMults); hold off;
xticks(linspace(1,24*pointsMultiplier,7));
xticklabels({'12:00 am','04:00 am','08:00 am','12:00 pm','04:00 pm','08:00 pm','12:00 am'});
title('All mutliplier values plotted'); legend('Residential Multipliers', 'Industrial Multipliers')
    
resDailyKwLoad = residentialMuls'*residentialKwLoad; % getting all 24*nptMult points for each houseLoad bus, kw
indDailyKwLoad = industrialMults'*industrialKwFromTestCase; % getting all 24*nptMult points for each industLoad bus, kw

resDailyKvarLoad = residentialMuls'*residentialKvarLoad; % getting all 24*nptMult points for each houseLoad bus, kvar
indDailyKvarLoad = industrialMults'*industrialKvarFromTestCase; % getting all 24*nptMult points for each industLoad bus, kvar

resRatedKwLoad = max(resDailyKwLoad, [], 1);
indRatedKwLoad = max(indDailyKwLoad, [], 1);

resRatedKvarLoad = max(resDailyKvarLoad, [], 1);
indRatedKvarLoad = max(indDailyKvarLoad, [], 1);

makingLoads(loadFileComponents, [industrialKwFromTestCase residentialKwLoad],...
    [industrialKvarFromTestCase residentialKvarLoad]);

% TODO: Change them to more meaningful names
resLoadsKwFinal = [];
resLoadsKvarFinal = [];
indLoadsKwFinal = [];
indLoadsKvar = [];

%% You can check the modified IEEE 13 test case in our thesis to verify bus names

% 20% of max load we'll set to rated; We hardcoded it -> Formula was this I suppose: ceil(1.2*max(sum(reskwDaily, 2) + sum(indkwDaily, 2)))
ratedSubStationKw = 4788; 

% All below are distribution trafo's in our model
distTrafokVA = 1000;
ratedDisTrafo634Kw = distTrafokVA * homePowerFactor; % 1000 is the distr. trafo kVA
ratedDisTrafo652Kw = distTrafokVA * homePowerFactor;
ratedDisTrafo611Kw = distTrafokVA * homePowerFactor;
ratedDisTrafo675Kw = distTrafokVA * homePowerFactor;
ratedDisTrafo646Kw = distTrafokVA * homePowerFactor;

%% Normal condition all plots

% Arrays that'll hold per-phase voltage values
% These will be cleared for each case
V1pu = []; % Phase A voltage in per unit (p.u.)
V2pu = []; % Phase B voltage in per unit (p.u.)
V3pu = []; % Phase C voltage in per unit (p.u.)

% power flows for substation and distribution transformer
% These will be cleared for each case
substationDailyKw = [];
substationDailyvar = [];

distTrafo652DailyKw = [];
distTrafo652DailyKvar = [];

distTrafo634DailyKw = [];
distTrafo634DailyKvar = [];

distTrafo611DailyKw = [];
distTrafo611DailyKvar = [];

distTrafo675DailyKw = [];
distTrafo675DailyKvar = [];

distTrafo646DailyKw = [];
distTrafo646DailyKvar = [];

substationCurrent = [];
distTrafo634LVCurrent = [];
distTrafo675LVCurrent = [];
distTrafo652LVCurrent = [];
distTrafo611LVCurrent = [];
distTrafo646LVCurrent = [];
bus680OutgoingCurrent = [];
bus692OutgoingCurrent = [];

PLoss = []; 
QLoss = [];

% Declared once
showLossFor = {'Transformer.SUB';
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
plotComponentNames = ["Transformer.SUB";'Transformer.XFM1';'Transformer.XFM2';'Transformer.XFM3';'Transformer.XFM4';
'Transformer.XFM5';'Line.650632';'Line.632671';'Line.671680';'Line.632633';'Line.632645';'Line.671684';
'Line.671692'; 'Capacitor.CAP1'; 'Capacitor.CAP2'; 'Load.645A'; 'Load.645B'; 'Load.645C'; 'Load.671';
'Load.692'; 'Load.680A'; 'Load.680B'; 'Load.680C'; 'Load.634A'; 'Load.634B'; 'Load.634C';
'Load.646A'; 'Load.646B'; 'Load.652'; 'Load.611'; 'Load.675A'; 'Load.675B';'Load.675C'];

% This will be defined once as well
AllLosses = [];

startLoadAt = 16;
phaseWiseArray = [4 5 6 0 0 1 2 3 1 2 3 1 2 1 3 1 2 3];
nTrafos = 6; % 5 distr. trafo's and one substation trafo

% clearvars reskwDaily reskvarDaily indkwDaily indkvarDaily

for i = 1:nPoints
    %residential house load
    kwLoopRes = resDailyKwLoad(i,:); kvarLoopRes = resDailyKvarLoad(i,:); % The particular power values for the i-th hour
    %industrial load                                              
    kwLoopInd = indDailyKwLoad(i,:); kvarLoopInd = indDailyKvarLoad(i,:); % The particular power values for the i-th hour
    
    makingLoads(loadFileComponents, [kwLoopInd kwLoopRes],...
        [kvarLoopInd kvarLoopRes]); % make the loadss.txt file
    DailyLoadFlow; % compute the daily load flow from the files in OpenDSS
    V1pu = [V1pu; Vpu1hr(:,1)']; % Get the per-phase voltage values for phase a
    V2pu = [V2pu; Vpu1hr(:,2)']; % Get the per-phase voltage values for phase b
    V3pu = [V3pu; Vpu1hr(:,3)']; % Get the per-phase voltage values for phase c
    
    % TODO: Replace xlsread with read csv

    % Read the powers per phase from the file
    [tempTrafoKw, tempTrafoKvar, tempLoadKw, tempLoadKvar]...
         = getPs(xlsread('MasterIEEE13_EXP_P_BYPHASE.CSV',1,'E2:J34'),...
         startLoadAt, phaseWiseArray, nTrafos); % Extract the power values for all loads and trafos
    
    indLoadsKwFinal = [indLoadsKwFinal; sum(tempLoadKw(1:nindustrialLoads,:), 1)];
    indLoadsKvar = [indLoadsKvar; sum(tempLoadKvar(1:nindustrialLoads,:), 1)];
    
    resLoadsKwFinal = [resLoadsKwFinal; sum(tempLoadKw(nindustrialLoads+1:end,:), 1)];
    resLoadsKvarFinal = [resLoadsKvarFinal; sum(tempLoadKvar(nindustrialLoads+1:end,:), 1)];
    
    tempI = getCurr3ph1hr(xlsread('MasterIEEE13_EXP_CURRENTS.CSV',1,'L2:Q8'), 2, 7);
    % tempBus680I = getCurr1ph1hr(xlsread('MasterIEEE13_EXP_CURRENTS.CSV',1,'B23:E25'));
    tempBus692I = getCurr3ph1hr(xlsread('MasterIEEE13_EXP_CURRENTS.CSV',1,'L2:Q15'), 10 ,14);
    
    substationCurrent = [substationCurrent; tempI(1,:)]; % Need this for final
    distTrafo634LVCurrent = [distTrafo634LVCurrent; tempI(2,:)]; % Need this for final
    distTrafo675LVCurrent = [distTrafo675LVCurrent; tempI(3,:)]; % Need this for final
    distTrafo652LVCurrent = [distTrafo652LVCurrent; tempI(4,:)];
    distTrafo611LVCurrent = [distTrafo611LVCurrent; tempI(5,:)];
    distTrafo646LVCurrent = [distTrafo646LVCurrent; tempI(6,:)];
    bus680OutgoingCurrent = [bus680OutgoingCurrent; tempBus692I(1,:)]; % Need this for final
    bus692OutgoingCurrent = [bus692OutgoingCurrent; tempBus692I(end,:)]; % Need this for final
    

    % substation and distr. trafo.
    substationDailyKw = [substationDailyKw; tempTrafoKw(1,:)]; % Store the kW power flow for substation
    subkvar = [subkvar; tempTrafoKvar(1,:)]; % Store the kVAR power flow for substation
    
    distTrafo634DailyKw = [distTrafo634DailyKw; tempTrafoKw(2,:)]; % Store the kW power flow for dist. trafo. of 634
    distTrafo634DailyKvar = [distTrafo634DailyKvar; tempTrafoKvar(2,:)]; % Store the kVAR power flow for dist. trafo. of 634
    
    distTrafo675DailyKw = [distTrafo675DailyKw; tempTrafoKw(3,:)]; % Store the kW power flow for dist. trafo. of 675
    distTrafo675DailyKvar = [distTrafo675DailyKvar; tempTrafoKvar(3,:)]; % Store the kVAR power flow for dist. trafo. of 675
    
    distTrafo652DailyKw = [distTrafo652DailyKw; tempTrafoKw(4,:)]; % Store the kW power flow for dist. trafo. of 652
    distTrafo652DailyKvar = [distTrafo652DailyKvar; tempTrafoKvar(4,:)]; % Store the kVAR power flow for dist. trafo. of 652
        
    distTrafo611DailyKw = [distTrafo611DailyKw; tempTrafoKw(5,:)]; % Store the kW power flow for dist. trafo. of 611
    distTrafo611DailyKvar = [distTrafo611DailyKvar; tempTrafoKvar(5,:)]; % Store the kVAR power flow for dist. trafo. of 611
    
    distTrafo646DailyKw = [distTrafo646DailyKw; tempTrafoKw(6,:)]; % Store the kW power flow for dist. trafo. of 646
    distTrafo646DailyKvar = [distTrafo646DailyKvar; tempTrafoKvar(6,:)]; % Store the kVAR power flow for dist. trafo. of 646
    
    [PLossTemp, QLossTemp] = getLinePower(readtable("MasterIEEE13_EXP_LOSSES.CSV")); % Getting all loses for the
                                                                                     % i-th HOUR
    PLoss = [PLoss PLossTemp];
    QLoss = [QLoss QLossTemp];

end

save('normal','V1pu','V2pu','V3pu',...
    'distTrafo634DailyKw', 'distTrafo634DailyKvar', 'distTrafo675DailyKw', 'distTrafo675DailyKvar', 'distTrafo652DailyKw', 'distTrafo652DailyKvar',...
    'distTrafo611DailyKw', 'distTrafo611DailyKvar', 'distTrafo646DailyKw', 'distTrafo646DailyKvar',...
    'substationDailyKw','subkvar', 'resLoadsKwFinal','indLoadsKwFinal',...
    'resLoadsKvarFinal','indLoadsKvar','showLossFor',...
    'resRatedKwLoad','resRatedKvarLoad','indRatedKwLoad','indRatedKvarLoad',...
    'substationCurrent','distTrafo634LVCurrent','distTrafo646LVCurrent','distTrafo675LVCurrent','distTrafo652LVCurrent','distTrafo611LVCurrent','bus680OutgoingCurrent','bus692OutgoingCurrent',...
    'PLoss', 'QLoss');

clc;

%%
close all;

V1pu(V1pu==0) = NaN;
V2pu(V2pu==0) = NaN;
V3pu(V3pu==0) = NaN;

%%
% Total load (residential + industrial) graphs
figure('Renderer', 'painters', 'Position', [200 -100 1100 700]);
hold on; 
plot(sum(resLoadsKwFinal, 2) + sum(indLoadsKwFinal, 2),'Color','blue','Marker','o','MarkerFaceColor','white',...
    'linewidth',plotWidth);
hold off;

xticks(linspace(1,24*pointsMultiplier,7));
xticklabels({'12:00 am','04:00 am','08:00 am','12:00 pm','04:00 pm','08:00 pm','12:00 am'});
title('Daily load profile for all loads at normal conditions', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
xlim([1 24*pointsMultiplier]); grid on;
ylim([(floor(min(sum(resLoadsKwFinal, 2) + sum(indLoadsKwFinal, 2))/500)*500 - 500),...
    (ceil(max(sum(resLoadsKwFinal, 2) + sum(indLoadsKwFinal, 2))/500)*500)])

busVSet = [2 4 5 6 8 14]; % Only the bus voltages I need

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
    plot(1.05*ones(1,nPoints),'linewidth',plotWidth,'linestyle','--','color','black');
    plot(0.95*ones(1,nPoints),'linewidth',plotWidth,'linestyle','--','color','black');
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
    lgnd.FontSize = legendFontSize; clear str;
    grid on;
    
    title(strcat("Voltage variation for bus", " ", busNames(i),", ","at normal conditions (base kV_{\phi} ="," ",...
        num2str(round(perUnitBasekVs(i)/sqrt(3),3))," ","kV_{rms})")); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Voltage in PU', 'FontSize', axisFontSize);
    ax = gca; ax.FontSize = axisFontSize;
    
    xticks(linspace(1,24*pointsMultiplier,7));
    xticklabels({'12:00 am','04:00 am','08:00 am','12:00 pm','04:00 pm','08:00 pm','12:00 am'});
    ylim([min([V1pu(:,i)' V2pu(:,i)' V3pu(:,i)' 0.95])-0.01 1.1])
    xlim([1 pointsMultiplier*24]);
end

%%
% Plotting 24-hr flow at substation
figure('Renderer', 'painters', 'Position', [200 -100 1100 700])
plot(substationDailyKw(:,1),'linewidth',plotWidth); hold on;
plot(substationDailyKw(:,2),'linewidth',plotWidth+0.2); 
plot(substationDailyKw(:,3),'linewidth',plotWidth+0.4);
plot(sum(substationDailyKw, 2),'linewidth',plotWidth+0.6);
plot(ones(1,nPoints)*ratedSubStationKw,'linewidth',plotWidth);
hold off;

xticks(linspace(1,24*pointsMultiplier,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Substation (bus 650) active power flow at normal conditions', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Sum of" + newline + "all phases";
lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = legendFontSize; clearvars str;
ylim([floor(min(substationDailyKw, [], 'all'))/500/1.25 ceil(max([sum(substationDailyKw, 2); ratedSubStationKw]))*1.225])
xlim([1 24*pointsMultiplier]);
%% I think this plot was not needed for final eval
% % Plotting 24-hr flow at substation (kVA plot) -> use sqrt(-1) instead of
                                                    % i because I used it as a counter
% figure('Renderer', 'painters', 'Position', [200 -100 1100 700])
% plot(abs(subkw(:,1)+sqrt(-1)*subkvar(:,1)),'linewidth',plotWidth); hold on;
% plot(abs(subkw(:,1)+sqrt(-1)*subkvar(:,1)),'linewidth',plotWidth+0.2); 
% plot(abs(subkw(:,1)+sqrt(-1)*subkvar(:,1)),'linewidth',plotWidth+0.4);
% plot(abs(sum(subkw, 2)+sqrt(-1)*sum(subkvar, 2)),'linewidth',plotWidth+0.6);
% plot(ones(1,npts)*5500,'linewidth',plotWidth);
% hold off;
% 
% xticks(linspace(1,24*nptMult,7));
%     xticklabels({'12:00 am','04:00 am','08:00 am',...
%         '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
% grid on;
% title('Substation (bus 650) active power flow at normal conditions', 'FontSize', titleFontSize); xlabel('Time',...
%     'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
% ax = gca; ax.FontSize = axisFontSize;
% str = "Sum of" + newline + "all phases";
% lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
% lgnd.FontSize = lgndFontSize; clearvars str;
% % ylim([floor(min(subkw, [], 'all'))/500 ceil(max(sum(subkw, 2))/500)*500*1.25])
% xlim([1 24*nptMult]);

%%
% Plotting 24-hr flow at dist. trafo. 634
figure('Renderer', 'painters', 'Position', [200 -100 1100 700])
plot(distTrafo634DailyKw(:,1),'linewidth',plotWidth); hold on;
plot(distTrafo634DailyKw(:,2),'linewidth',plotWidth+0.2); 
plot(distTrafo634DailyKw(:,3),'linewidth',plotWidth+0.4); 
plot(sum(distTrafo634DailyKw, 2),'linewidth',plotWidth+0.6);
plot(ones(1,nPoints)*ratedDisTrafo634Kw,'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*pointsMultiplier,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Distribution transformer (634) power flow at normal conditions', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Sum of" + newline + "all phases";
ylim([floor(min(min(distTrafo634DailyKw), [], 'all')/500) ceil(max([sum(distTrafo634DailyKw, 2); ratedDisTrafo634Kw]))*1.225])
xlim([1 24*pointsMultiplier]);
lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = legendFontSize; clearvars str;

%%
% Plotting 24-hr flow at dist. trafo. 652
figure('Renderer', 'painters', 'Position', [200 -100 1100 700])
plot(distTrafo652DailyKw(:,1),'linewidth',plotWidth*3); hold on;
plot(distTrafo652DailyKw(:,2),'linewidth',plotWidth*1.5); 
plot(distTrafo652DailyKw(:,3),'linewidth',plotWidth*0.8); 
plot(sum(distTrafo652DailyKw, 2),'linewidth',plotWidth*1.15);
plot(ones(1,nPoints)*ratedDisTrafo652Kw,'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*pointsMultiplier,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Distribution transformer (652) power flow at normal conditions', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Sum of" + newline + "all phases";
ylim([floor(min(min(distTrafo652DailyKw), [], 'all')/500) ceil(max([sum(distTrafo652DailyKw, 2); ratedDisTrafo652Kw]))*1.225])
xlim([1 24*pointsMultiplier]);
lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = legendFontSize; clearvars str;

%%
% Plotting 24-hr flow at dist. trafo. 611
figure('Renderer', 'painters', 'Position', [200 100 1100 700])
plot(distTrafo611DailyKw(:,1),'linewidth',plotWidth); hold on;
plot(distTrafo611DailyKw(:,2),'linewidth',plotWidth*1.12); 
plot(distTrafo611DailyKw(:,3),'linewidth',plotWidth*1.14); 
plot(sum(distTrafo611DailyKw, 2),'linewidth',plotWidth*1.16);
plot(ones(1,nPoints)*ratedDisTrafo611Kw,'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*pointsMultiplier,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Distribution transformer (611) power flow at normal conditions', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Sum of" + newline + "all phases";
ylim([floor(min(min(distTrafo611DailyKw), [], 'all')/500) ceil(max([sum(distTrafo611DailyKw, 2); ratedDisTrafo611Kw]))*1.225])
xlim([1 24*pointsMultiplier]);
lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = legendFontSize; clearvars str;

%%
% Plotting 24-hr flow at dist. trafo. 675
figure('Renderer', 'painters', 'Position', [200 100 1100 700])
plot(distTrafo675DailyKw(:,1),'linewidth',plotWidth); hold on;
plot(distTrafo675DailyKw(:,2),'linewidth',plotWidth+0.2); 
plot(distTrafo675DailyKw(:,3),'linewidth',plotWidth+0.4); 
plot(sum(distTrafo675DailyKw, 2),'linewidth',plotWidth+0.6);
plot(ones(1,nPoints)*ratedDisTrafo675Kw,'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*pointsMultiplier,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Distribution transformer (675) power flow at normal conditions', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Sum of" + newline + "all phases";
ylim([floor(min(min(distTrafo675DailyKw), [], 'all')/500) ceil(max([sum(distTrafo675DailyKw, 2); ratedDisTrafo652Kw]))*1.225])
xlim([1 24*pointsMultiplier]);
lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = legendFontSize; clearvars str;

%%
% Plotting 24-hr flow at dist. trafo. 646
figure('Renderer', 'painters', 'Position', [200 100 1100 700])
plot(distTrafo646DailyKw(:,1),'linewidth',plotWidth); hold on;
plot(distTrafo646DailyKw(:,2),'linewidth',plotWidth*1.12); 
plot(distTrafo646DailyKw(:,3),'linewidth',plotWidth*1.14); 
plot(sum(distTrafo646DailyKw, 2),'linewidth',plotWidth*1.16);
plot(ones(1,nPoints)*ratedDisTrafo646Kw,'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*pointsMultiplier,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Distribution transformer (646) power flow at normal conditions', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Sum of" + newline + "all phases";
ylim([floor(min(min(distTrafo646DailyKw), [], 'all')/500) ceil(max([sum(distTrafo646DailyKw, 2); ratedDisTrafo646Kw]))*1.225])
xlim([1 24*pointsMultiplier]);
lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = legendFontSize; clearvars str;

%%
% Plotting 24-hr current flow at 2ndary terminal for sub
figure('Renderer', 'painters', 'Position', [200 100 1100 700])
plot(substationCurrent(:,1),'linewidth',plotWidth); hold on;
plot(substationCurrent(:,2),'linewidth',plotWidth*1.12); 
plot(substationCurrent(:,3),'linewidth',plotWidth*1.14); 
% plot(sum(distTrafo646DailyKw, 2),'linewidth',plotWidth*1.16);
plot(ones(1,nPoints)*(1.1*5500/sqrt(3)/11),'linewidth',plotWidth);
% Rated value for "Dove" conductor should not be used with 
% a trafo
% plot(ones(1,npts)*(5500*1.5/11/sqrt(3)),'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*pointsMultiplier,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Substation (Sub) current flow at normal conditions', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Current (A)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Normal Ampere rating";
ylim([floor(min((substationCurrent), [], 'all'))/500 ceil(max([max(substationCurrent,[] ,'all'); 1.1*5500/sqrt(3)/11]))*1.15]) % Max rating are from max current,
                                                                                                                       % normamps and emergamps
xlim([1 24*pointsMultiplier]);
lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = legendFontSize; clearvars str;

%%
% Plotting 24-hr current flow at 2ndary terminal for 634
figure('Renderer', 'painters', 'Position', [200 100 1100 700])
plot(distTrafo634LVCurrent(:,1),'linewidth',plotWidth); hold on;
plot(distTrafo634LVCurrent(:,2),'linewidth',plotWidth*1.12); 
plot(distTrafo634LVCurrent(:,3),'linewidth',plotWidth*1.14); 
% plot(sum(distTrafo646DailyKw, 2),'linewidth',plotWidth*1.16);
plot(ones(1,nPoints)*1000*1.1/0.4/sqrt(3),'linewidth',plotWidth);
% plot(ones(1,npts)*(1000*1.5/0.4/sqrt(3)),'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*pointsMultiplier,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Distribution transformer (634) current flow at normal conditions', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Current (A)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Normal Ampere rating";
ylim([floor(min((distTrafo634LVCurrent), [], 'all'))/500 ceil(max([max(distTrafo634LVCurrent,[] ,'all'); 1000*1.1/0.4/sqrt(3)]))*1.15]) % Max rating are from max current,
                                                                                                      % normamps 
xlim([1 24*pointsMultiplier]);
lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = legendFontSize; clearvars str;

%%
% Plotting 24-hr current flow at 2ndary terminal for 675
figure('Renderer', 'painters', 'Position', [200 100 1100 700])
plot(distTrafo675LVCurrent(:,1),'linewidth',plotWidth); hold on;
plot(distTrafo675LVCurrent(:,2),'linewidth',plotWidth*1.12); 
plot(distTrafo675LVCurrent(:,3),'linewidth',plotWidth*1.14); 
% plot(sum(distTrafo646DailyKw, 2),'linewidth',plotWidth*1.16);
plot(ones(1,nPoints)*1000*1.1/0.4/sqrt(3),'linewidth',plotWidth);
% plot(ones(1,npts)*(1000*1.5/0.4/sqrt(3)),'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*pointsMultiplier,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Distribution transformer (675) current flow at normal conditions', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Current (A)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Normal Ampere rating";
ylim([floor(min((distTrafo675LVCurrent), [], 'all'))/500 ceil(max([max(distTrafo675LVCurrent,[] ,'all'); 1000*1.1/0.4/sqrt(3)]))*1.15]) % Max rating are from max current,
                                                                                                      % normamps 
xlim([1 24*pointsMultiplier]);
lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = legendFontSize; clearvars str;

%%
% Plotting 24-hr current flow at 2ndary terminal for bus 680
% IMPORTANT! I'll mention the current rating separately
figure('Renderer', 'painters', 'Position', [200 100 1100 700])
plot(bus680OutgoingCurrent(:,1),'linewidth',plotWidth); hold on;
plot(bus680OutgoingCurrent(:,2),'linewidth',plotWidth*1.12); 
plot(bus680OutgoingCurrent(:,3),'linewidth',plotWidth*1.14); 
% plot(ones(1,npts)*726,'linewidth',plotWidth); % CURRENT RATING!!!
% plot(ones(1,npts)*(1000*1.5/0.4/sqrt(3)),'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*pointsMultiplier,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Bus 680 current flow at normal conditions', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Current (A)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Normal Ampere rating";
ylim([floor(min((bus680OutgoingCurrent), [], 'all'))/500 ceil(max([max(bus680OutgoingCurrent,[] ,'all'); 0]))*1.15]) % Max rating are from max current,
                                                                                                      % normamps 
xlim([1 24*pointsMultiplier]);
lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = legendFontSize; clearvars str;

%%
% Plotting 24-hr current flow at 2ndary terminal for bus 692
% IMPORTANT! I'll mention the current rating separately
figure('Renderer', 'painters', 'Position', [200 100 1100 700])
plot(bus692OutgoingCurrent(:,1),'linewidth',plotWidth); hold on;
plot(bus692OutgoingCurrent(:,2),'linewidth',plotWidth*1.12); 
plot(bus692OutgoingCurrent(:,3),'linewidth',plotWidth*1.14); 
% plot(ones(1,npts)*260,'linewidth',plotWidth); % CURRENT RATING!!!
% We'll just mention it!!!
% plot(ones(1,npts)*(1000*1.5/0.4/sqrt(3)),'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*pointsMultiplier,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Bus 692 current flow at normal conditions', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Current (A)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Normal Ampere rating";
ylim([floor(min((bus692OutgoingCurrent), [], 'all'))/500 ceil(max([max(bus692OutgoingCurrent,[] ,'all'); 0]))*1.15]) % Max rating are from max current,
                                                                                                      % normamps 
xlim([1 24*pointsMultiplier]);
lgnd = legend({'Phase a','Phase b','Phase c'},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = legendFontSize; clearvars str;

%%
% Plotting 24-hr current flow at 2ndary terminal for bus 652
figure('Renderer', 'painters', 'Position', [200 100 1100 700])
plot(distTrafo652LVCurrent(:,1),'linewidth',plotWidth); hold on;
plot(distTrafo652LVCurrent(:,2),'linewidth',plotWidth*1.12); 
plot(distTrafo652LVCurrent(:,3),'linewidth',plotWidth*1.14); 
plot(ones(1,nPoints)*1000*1.1/0.4/sqrt(3),'linewidth',plotWidth);
% We'll just mention it!!!
% plot(ones(1,npts)*(1000*1.5/0.4/sqrt(3)),'linewidth',plotWidth);
hold off;
xticks(linspace(1,24*pointsMultiplier,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title('Distribution transformer (652) current flow at normal conditions', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Current (A)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str = "Normal Ampere rating";
ylim([floor(min((distTrafo652LVCurrent), [], 'all'))/500 ceil(max([max(distTrafo652LVCurrent,[] ,'all'); 0]))*1.15]) % Max rating are from max current,
                                                                                                      % normamps 
xlim([1 24*pointsMultiplier]);
lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southoutside','orientation','horizontal');
lgnd.FontSize = legendFontSize; clearvars str;

%% I think following plots were not needed for final eval
% % Plotting 24-hr current flow at 2ndary terminal for 675
% figure('Renderer', 'painters', 'Position', [200 100 1100 700])
% plot(distTrafo675LVCurrent(:,1),'linewidth',plotWidth); hold on;
% plot(distTrafo675LVCurrent(:,2),'linewidth',plotWidth*1.12); 
% plot(distTrafo675LVCurrent(:,3),'linewidth',plotWidth*1.14); 
% % plot(sum(distTrafo646DailyKw, 2),'linewidth',plotWidth*1.16);
% plot(ones(1,npts)*(357),'linewidth',plotWidth);
% % plot(ones(1,npts)*(1000*1.5/0.4/sqrt(3)),'linewidth',plotWidth);
% hold off;
% xticks(linspace(1,24*nptMult,7));
%     xticklabels({'12:00 am','04:00 am','08:00 am',...
%         '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
% grid on;
% title('Distribution transformer (634) current flow at normal conditions', 'FontSize', titleFontSize); xlabel('Time',...
%     'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
% ax = gca; ax.FontSize = axisFontSize;
% str = "Normal Ampere flow";
% ylim([floor(min((distTrafo675LVCurrent), [], 'all'))/500 ceil(max([sum(distTrafo646DailyKw, 2); 357]))*1.15]) % Max rating are from max current,
%                                                                                   % normamps and emergamps
% xlim([1 24*nptMult]);
% lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southoutside','orientation','horizontal');
% lgnd.FontSize = lgndFontSize; clearvars str;
% 
% %%
% % Plotting 24-hr current flow at 2ndary terminal for 634
% figure('Renderer', 'painters', 'Position', [200 100 1100 700])
% plot(bus680OutgoingCurrent(:,1),'linewidth',plotWidth); hold on;
% plot(bus680OutgoingCurrent(:,2),'linewidth',plotWidth*1.12); 
% plot(bus680OutgoingCurrent(:,3),'linewidth',plotWidth*1.14); 
% % plot(sum(distTrafo646DailyKw, 2),'linewidth',plotWidth*1.16);
% % plot(ones(1,npts)*(5500*1.1/11/sqrt(3)),'linewidth',plotWidth);
% % plot(ones(1,npts)*(5500*1.5/11/sqrt(3)),'linewidth',plotWidth);
% % We need to find SOME rated current value and I think I can find it from
% % OpenDSS
% hold off;
% xticks(linspace(1,24*nptMult,7));
%     xticklabels({'12:00 am','04:00 am','08:00 am',...
%         '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
% grid on;
% title('Distribution transformer (Sub) current flow at normal conditions', 'FontSize', titleFontSize); xlabel('Time',...
%     'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
% ax = gca; ax.FontSize = axisFontSize;
% str = "Normal Ampere flow";
% ylim([floor(min((substationCurrent), [], 'all'))/500 ceil(max([sum(distTrafo646DailyKw, 2); 5500*1.1/11/sqrt(3); 5500*1.5/11/sqrt(3)]))*1.15]) % Max rating are from max current,
%                                                                                                                        % normamps and emergamps
% xlim([1 24*nptMult]);
% lgnd = legend({'Phase a','Phase b','Phase c',str,"Emergency Ampere flow"},'Location','southoutside','orientation','horizontal');
% lgnd.FontSize = lgndFontSize; clearvars str;

