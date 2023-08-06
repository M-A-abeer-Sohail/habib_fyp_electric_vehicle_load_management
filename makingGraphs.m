%% Set up
clear;
close all;

nptMult = 4;
npts = 24*nptMult;

color1 = [34;139;34]'/255;

ratedSubkw = 4788; % 20% of max load we'll set to rated.
ratedDisTr634kw = 617;
ratedDisTr652kw = 380;
ratedDisTr611kw = 374;
ratedDisTr675kw = 1086;
ratedDisTr646kw = 552;

axisFontSize = 22;
lgndFontSize = 14;
plotWidth = 5;
titleFontSize = 16;

normalCond = load("C:\Users\Abeer\OneDrive - Habib University\V2G Potential FYP\FYP codes & data\Capstone II\03e running the EV loads\normal.mat");
ev40Uncoord = load("C:\Users\Abeer\OneDrive - Habib University\V2G Potential FYP\FYP codes & data\Capstone II\03e running the EV loads\ev40Pent.mat");
ev90Uncoord = load("C:\Users\Abeer\OneDrive - Habib University\V2G Potential FYP\FYP codes & data\Capstone II\03e running the EV loads\ev90Pent.mat");
ev40Coord = load('ev40Pent.mat');
ev90Coord = load('ev90Pent.mat');

load busNames.mat;
basekVs = [132 11 11 0.4 11 0.4 11 0.4 0.4 11 0.4 11 11 11];

%% Plotting daily load flow

% [178 153 0] / 255% the color [1.0000    0.4392         0]  [0    0.4863    1.0000]

axisFontSize = 22;
lgndFontSize = 12;
plotWidth = 5;
titleFontSize = 16;

close all;
% Total load (residential + industrial + ev) graphs
figure('Renderer', 'painters', 'Position', [200 -100 1100 700]);
t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile;
hold on;
plot(sum(normalCond.resLdkw, 2) + sum(normalCond.industLdkw, 2),'color',[0    0.6510    0.1569],...
    'linewidth',plotWidth-1);
plot(sum(ev40Uncoord.PevTot, 2)+sum(ev40Uncoord.resLdkw, 2) + sum(ev40Uncoord.industLdkw, 2),'Color','red','linestyle','-.',...
    'linewidth',plotWidth);
plot(sum(ev90Uncoord.PevTot, 2)+ sum(ev90Uncoord.resLdkw, 2) + sum(ev90Uncoord.industLdkw, 2),'Color','blue','linestyle','-.',...
    'linewidth',plotWidth);
plot(sum(ev40Coord.PevTot, 2)+sum(ev40Coord.resLdkw, 2) + sum(ev40Coord.industLdkw, 2),'Color','red',...
    'linewidth',plotWidth);
plot(sum(ev90Coord.PevTot, 2)+sum(ev90Coord.resLdkw, 2) + sum(ev90Coord.industLdkw, 2),'Color','blue',...
    'linewidth',plotWidth);
hold off;

xticks(linspace(1,24*nptMult,7));
xticklabels({'12:00 am','04:00 am','08:00 am','12:00 pm','04:00 pm','08:00 pm','12:00 am'});
title('Daily load profile of uncoordinated vs. coordinated schemes', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
xlim([1 24*nptMult]); grid on;
ylim([(floor(min([sum(ev40Uncoord.PevTot, 2)+sum(ev40Uncoord.resLdkw, 2) + sum(ev40Uncoord.industLdkw, 2)...
    sum(ev90Uncoord.PevTot, 2)+ sum(ev90Uncoord.resLdkw, 2) + sum(ev90Uncoord.industLdkw, 2)...
    sum(ev40Coord.PevTot, 2)+sum(ev40Coord.resLdkw, 2) + sum(ev40Coord.industLdkw, 2)...
    sum(ev90Coord.PevTot, 2)+sum(ev90Coord.resLdkw, 2) + sum(ev90Coord.industLdkw, 2)],[],'all')/500)*500*0.95),...
    (ceil(max([sum(ev40Uncoord.PevTot, 2)+sum(ev40Uncoord.resLdkw, 2) + sum(ev40Uncoord.industLdkw, 2)...
    sum(ev90Uncoord.PevTot, 2)+ sum(ev90Uncoord.resLdkw, 2) + sum(ev90Uncoord.industLdkw, 2)...
    sum(ev40Coord.PevTot, 2)+sum(ev40Coord.resLdkw, 2) + sum(ev40Coord.industLdkw, 2)...
    sum(ev90Coord.PevTot, 2)+sum(ev90Coord.resLdkw, 2) + sum(ev90Coord.industLdkw, 2)],[],'all')/500)*500*1.0625)])
tempstr = ["Normal conditions",string(strjoin(["40% EV" newline "penetration with" newline "uncoordinated charging"])),...
           string(strjoin(["90% EV" newline "penetration with" newline "uncoordinated charging"])),...
           string(strjoin(["40% EV" newline "penetration with" newline "coordinated charging"])),...
           string(strjoin(["90% EV" newline "penetration with" newline "coordinated charging"]))];
lgnd = legend({tempstr(1), tempstr(2), tempstr(3), tempstr(4), tempstr(5)},...
    'Location','northwest','orientation','vertical');
    lgnd.FontSize = lgndFontSize;
    
%% Plotting total Non EV loads
close all;
% Total load (residential + industrial + ev) graphs
figure('Renderer', 'painters', 'Position', [200 -100 1100 700]);
t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile;
hold on; 
plot(sum(ev40Uncoord.resLdkw, 2) + sum(ev40Uncoord.industLdkw, 2),'Color','blue',...
    'linewidth',plotWidth);
plot(sum(ev90Uncoord.resLdkw, 2) + sum(ev90Uncoord.industLdkw, 2),'Color',color1,...
    'linewidth',plotWidth);
plot(sum(ev40Coord.resLdkw, 2) + sum(ev40Coord.industLdkw, 2),'Color',[0.4940, 0.1840, 0.5560],...
    'linewidth',plotWidth);
plot(sum(ev90Coord.resLdkw, 2) + sum(ev90Coord.industLdkw, 2),'Color','red',...
    'linewidth',plotWidth);
hold off;

xticks(linspace(1,24*nptMult,7));
xticklabels({'12:00 am','04:00 am','08:00 am','12:00 pm','04:00 pm','08:00 pm','12:00 am'});
title('Daily load profile for all Non EV loads in uncoordinated charging', 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
xlim([1 24*nptMult]); grid on;
ylim([(floor(min(sum(ev40Uncoord.resLdkw, 2) + sum(ev40Uncoord.industLdkw, 2))/500)*500*0.95),...
    (ceil(max(sum(ev90Coord.PevTot, 2)+sum(ev40Uncoord.resLdkw, 2) + sum(ev40Uncoord.industLdkw, 2))/500)*500*0.95)])
tempstr = string(strjoin(["Normal conditions" newline "without EVs"]));
lgnd = legend({tempstr, '20% EV penetration', '40% EV penetration', '90% EV penetration'},...
    'Location','southoutside','orientation','horizontal');
    lgnd.FontSize = lgndFontSize;

%% Plotting the bus voltages of each bus 2x2
close all;
busVSet = [4]; % Only the bus voltages I need
axisFontSize = 22;
lgndFontSize = 15;
plotWidth = 5;
titleFontSize = 2;
sgTitleFontSize = 26;
close all;

strTrafo = ["3 phase primary substation transformer", "3 phase distribution transformer at bus 675", "3 phase distribution transformer at bus 652"];

for j = 2:length(busNames) % ignoring SOURCEBUS
    i = j;% busVSet(j);
    figure('Renderer', 'painters', 'Position', [200 -100 1300 900])
    t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
%     sgtitle("Voltage variation for " + strTrafo(j) + " (base kV = " +  num2str(basekVs(i)) + " kV)",'fontweight','bold',...
%         'fontsize',sgTitleFontSize);

    sgtitle("Voltage variation for " + busNames(i) + " (base kV = " +  num2str(basekVs(i)) + " kV)" + newline + "in uncoordinated vs coordinated schemes",'fontweight','bold',...
         'fontsize',sgTitleFontSize);
     
    legndArr = {};
    if ~isnan(ev90Uncoord.V1pu(1,i))
        legndArr(end+1) = {"Phase a"};
    end
    if ~isnan(ev90Uncoord.V2pu(1,i))
         legndArr(end+1) = {"Phase b"};
    end
    if ~isnan(ev90Uncoord.V3pu(1,i))
         legndArr(end+1) = {"Phase c"};
    end
    
    str1 = "1.05 p.u." + newline + "limit";
    str2 = "0.95 p.u." + newline + "limit";
    % legndArr(end+1) = {str1};
    legndArr(end+1) = {str2};
    clearvars str1 str2
    if ~isnan(normalCond.V1pu(1,i))
        legndArr(end+1) = {"Phase a in" + newline + "normal conditions"};
    end
    if ~isnan(normalCond.V2pu(1,i))
         legndArr(end+1) = {"Phase b in" + newline + "normal conditions"};
    end
    if ~isnan(normalCond.V3pu(1,i))
         legndArr(end+1) = {"Phase c in" + newline + "normal conditions"};
    end
    
    % 40% penetration uncoord
    nexttile;
    hold on;
    if ~isnan(ev40Uncoord.V1pu(1,i))
         plot(ev40Uncoord.V1pu(:,i),'linewidth',plotWidth);
    end
    if ~isnan(ev40Uncoord.V2pu(1,i))
         plot(ev40Uncoord.V2pu(:,i),'linewidth',plotWidth);
    end
    if ~isnan(ev40Uncoord.V3pu(1,i))
         plot(ev40Uncoord.V3pu(:,i),'linewidth',plotWidth);
    end
    % plot(1.05*ones(1,npts),'linewidth',plotWidth,'linestyle','--','color','black');
    plot(0.95*ones(1,npts),'linewidth',plotWidth,'linestyle','--','color','black');
    grid on;
   
    if ~isnan(normalCond.V1pu(1,i))
         plot(normalCond.V1pu(:,i),'linewidth',plotWidth-1.5,'color',[0, 0.5, 0],...
             'linestyle','-.');
    end
    if ~isnan(normalCond.V2pu(1,i))
         plot(normalCond.V2pu(:,i),'linewidth',plotWidth-1.5,'color',[0.75, 0, 0.75],...
             'linestyle','-.');
    end
    if ~isnan(normalCond.V3pu(1,i))
         plot(normalCond.V3pu(:,i),'linewidth',plotWidth-1.5,'color',[0.6350, 0.0780, 0.1840],...
             'linestyle','-.');
    end
   
    hold off; 
    
    titleStr = "40% EV penetration" + newline + "with uncoordinated charging";
    
    title((titleStr),'fontsize',titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Voltage in PU', 'FontSize', axisFontSize);
    ax = gca; ax.FontSize = axisFontSize;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am','12:00 pm','12:00 am'});
    ylim([min([ev40Uncoord.V1pu(:,i)' ev40Uncoord.V2pu(:,i)' ev40Uncoord.V3pu(:,i)' 0.95])-0.005 max([normalCond.V1pu(:,i)' normalCond.V2pu(:,i)' normalCond.V3pu(:,i)' 1])])
    xlim([1 nptMult*24]);
    
    % 90% Penetration uncoord
    nexttile;
    hold on;
    if ~isnan(ev40Coord.V1pu(1,i))
         plot(ev40Coord.V1pu(:,i),'linewidth',plotWidth);
    end
    if ~isnan(ev40Coord.V2pu(1,i))
         plot(ev40Coord.V2pu(:,i),'linewidth',plotWidth);
    end
    if ~isnan(ev40Coord.V3pu(1,i))
         plot(ev40Coord.V3pu(:,i),'linewidth',plotWidth);
    end
    % plot(1.05*ones(1,npts),'linewidth',plotWidth,'linestyle','--','color','black');
    plot(0.95*ones(1,npts),'linewidth',plotWidth,'linestyle','--','color','black');
    
    if ~isnan(normalCond.V1pu(1,i))
         plot(normalCond.V1pu(:,i),'linewidth',plotWidth-1.5,'color',[0, 0.5, 0],...
             'linestyle','-.');
    end
    if ~isnan(normalCond.V2pu(1,i))
         plot(normalCond.V2pu(:,i),'linewidth',plotWidth-1.5,'color',[0.75, 0, 0.75],...
             'linestyle','-.');
    end
    if ~isnan(normalCond.V3pu(1,i))
         plot(normalCond.V3pu(:,i),'linewidth',plotWidth-1.5,'color',[0.6350, 0.0780, 0.1840],...
             'linestyle','-.');
    end
    
    hold off;
    
    grid on;
    
        titleStr = "40% EV penetration" + newline + "with coordinated charging";
    
    title((titleStr),'fontsize',titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Voltage in PU', 'FontSize', axisFontSize);
    ax = gca; ax.FontSize = axisFontSize;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am','12:00 pm','12:00 am'});
    ylim([min([ev40Coord.V1pu(:,i)' ev40Coord.V2pu(:,i)' ev40Coord.V3pu(:,i)' 0.95])-0.005 max([normalCond.V1pu(:,i)' normalCond.V2pu(:,i)' normalCond.V3pu(:,i)' 1])])
    xlim([1 nptMult*24]);
    
    % 40% Penetration uncoord
    nexttile;
    hold on;
    if ~isnan(ev90Uncoord.V1pu(1,i))
         plot(ev90Uncoord.V1pu(:,i),'linewidth',plotWidth);
    end
    if ~isnan(ev90Uncoord.V2pu(1,i))
         plot(ev90Uncoord.V2pu(:,i),'linewidth',plotWidth);
    end
    if ~isnan(ev90Uncoord.V3pu(1,i))
         plot(ev90Uncoord.V3pu(:,i),'linewidth',plotWidth);
    end
    % plot(1.05*ones(1,npts),'linewidth',plotWidth,'linestyle','--','color','black');
    plot(0.95*ones(1,npts),'linewidth',plotWidth,'linestyle','--','color','black');
    
    if ~isnan(normalCond.V1pu(1,i))
         plot(normalCond.V1pu(:,i),'linewidth',plotWidth-1.5,'color',[0, 0.5, 0],...
             'linestyle','-.');
    end
    if ~isnan(normalCond.V2pu(1,i))
         plot(normalCond.V2pu(:,i),'linewidth',plotWidth-1.5,'color',[0.75, 0, 0.75],...
             'linestyle','-.');
    end
    if ~isnan(normalCond.V3pu(1,i))
         plot(normalCond.V3pu(:,i),'linewidth',plotWidth-1.5,'color',[0.6350, 0.0780, 0.1840],...
             'linestyle','-.');
    end
    
    hold off;

    grid on;
    
    titleStr = "90% EV penetration" + newline + "with uncoordinated charging";
    
    title((titleStr),'fontsize',titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Voltage in PU', 'FontSize', axisFontSize);
    ax = gca; ax.FontSize = axisFontSize;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am','12:00 pm','12:00 am'});
    ylim([min([ev90Uncoord.V1pu(:,i)' ev90Uncoord.V2pu(:,i)' ev90Uncoord.V3pu(:,i)' 0.95])-0.005 max([normalCond.V1pu(:,i)' normalCond.V2pu(:,i)' normalCond.V3pu(:,i)' 1])])
    xlim([1 nptMult*24]);
    
    grid on;
    
    % 90 % Penetration
    nexttile;
    hold on;
    if ~isnan(ev90Coord.V1pu(1,i))
         plot(ev90Coord.V1pu(:,i),'linewidth',plotWidth);
    end
    if ~isnan(ev90Coord.V2pu(1,i))
         plot(ev90Coord.V2pu(:,i),'linewidth',plotWidth);
    end
    if ~isnan(ev90Coord.V3pu(1,i))
         plot(ev90Coord.V3pu(:,i),'linewidth',plotWidth);
    end
    % plot(1.05*ones(1,npts),'linewidth',plotWidth,'linestyle','--','color','black');
    plot(0.95*ones(1,npts),'linewidth',plotWidth,'linestyle','--','color','black');
    
    if ~isnan(normalCond.V1pu(1,i))
         plot(normalCond.V1pu(:,i),'linewidth',plotWidth-1.5,'color',[0, 0.5, 0],...
             'linestyle','-.');
    end
    if ~isnan(normalCond.V2pu(1,i))
         plot(normalCond.V2pu(:,i),'linewidth',plotWidth-1.5,'color',[0.75, 0, 0.75],...
             'linestyle','-.');
    end
    if ~isnan(normalCond.V3pu(1,i))
         plot(normalCond.V3pu(:,i),'linewidth',plotWidth-1.5,'color',[0.6350, 0.0780, 0.1840],...
             'linestyle','-.');
    end
    
    hold off;
     
    grid on;
    titleStr = "90% EV penetration" + newline + "with coordinated charging";
    
    title((titleStr),'fontsize',titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Voltage in PU', 'FontSize', axisFontSize);
    ax = gca; ax.FontSize = axisFontSize;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am','12:00 pm','12:00 am'});
    ylim([min([ev90Coord.V1pu(:,i)' ev90Coord.V2pu(:,i)' ev90Coord.V3pu(:,i)' 0.95])-0.005 max([normalCond.V1pu(:,i)' normalCond.V2pu(:,i)' normalCond.V3pu(:,i)' 1])])
    xlim([1 nptMult*24]);
    
    % lgnd = legend(string(legndArr),'location','southeastoutside',...
    % 'Orientation','horizontal',...
    % 'FontSize',lgndFontSize);
end

%% Plot currents

% busNameSet = busVSet(1);
close all
axisFontSize = 22;
lgndFontSize = 14;
plotWidth = 5;
titleFontSize = 12;
sgTitleFontSize = 22;
close all;

sgtitleText = "Current variation for distribution transformer at bus 652 (base kV = 0.4 kV)";
plot1Val = ev40Uncoord.xfm652I;
plot2Val = ev40Coord.xfm652I;
plot3Val = ev90Uncoord.xfm652I;
plot4Val = ev90Coord.xfm652I;
commonVal = normalCond.xfm652I;

tileTitle1 = "40% EV penetration" + newline + "with uncoordinated charging";
tileTitle2 = "40% EV penetration" + newline + "with coordinated charging";
tileTitle3 = "90% EV penetration" + newline + "with uncoordinated charging";
tileTitle4 = "90% EV penetration" + newline + "with coordinated charging";

ratedVal = 1.1*1000/sqrt(3)/0.4;
ratedCol = [125, 92, 14]/255;

legndArr = {};

    if (plot1Val(1,1)>=1e-5)
        legndArr(end+1) = {"Phase a"};
    end
    if (plot1Val(1,2)>=1e-5)
         legndArr(end+1) = {"Phase b"};
    end
    if (plot1Val(1,3)>=1e-5)
         legndArr(end+1) = {"Phase c"};
    end
    
    str1 = "Normal Ampere" + newline + "rating";
    % legndArr(end+1) = {str1};
    legndArr(end+1) = {str1};
    clearvars str1 str2

    if (commonVal(1,1)>=1e-5)
        legndArr(end+1) = {"Phase a in" + newline + "normal conditions"};
    end
    if (commonVal(1,2)>=1e-5)
         legndArr(end+1) = {"Phase b in" + newline + "normal conditions"};
    end
    if (commonVal(1,3)>=1e-5)
         legndArr(end+1) = {"Phase c in" + newline + "normal conditions"};
    end
    
    figure('Renderer', 'painters', 'Position', [200 -100 1300 900])
    t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
    sgtitle(sgtitleText,'fontweight','bold',...
        'fontsize',sgTitleFontSize);
    
nexttile
hold on;
    if (plot1Val(1,1)>=1e-5)
         plot(plot1Val(:,1),'linewidth',plotWidth);
    end
    if (plot1Val(1,2)>=1e-5)
         plot(plot1Val(:,2),'linewidth',plotWidth);
    end
    if (plot1Val(1,3)>=1e-5)
         plot(plot1Val(:,3),'linewidth',plotWidth);
    end
    plot(ones(1,npts)*ratedVal,'linewidth',plotWidth,'color',ratedCol);
    grid on;
    
    if (commonVal(1,1)>=1e-5)
         plot(commonVal(:,1),'linewidth',plotWidth,'color',[0, 0.5, 0],...
             'linestyle','-.');
    end
    if (commonVal(1,2)>=1e-5)
         plot(commonVal(:,2),'linewidth',plotWidth,'color',[0.75, 0, 0.75],...
         'linestyle','-.');
    end
    if (commonVal(1,3)>=1e-5)
         plot(commonVal(:,3),'linewidth',plotWidth,'color',[0.6350, 0.0780, 0.1840],...
             'linestyle','-.');
    end
    hold off;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am','12:00 pm','12:00 am'});
    
    ylim([floor(min((commonVal), [], 'all'))*0.95 ceil(max([max(plot1Val,[] ,'all'); ratedVal]))*1.0375]) % Max rating are from max current,
                                                                                                                        % normamps and emergamps
    
    title(tileTitle1,'FontSize', titleFontSize);
    xlim([1 24*nptMult]);
    ax = gca; ax.FontSize = axisFontSize;
    
    xlabel('Time','FontSize', axisFontSize); ylabel('Current (A)','FontSize', axisFontSize);

nexttile
hold on;

    if (plot2Val(1,1)>=1e-5)
         plot(plot2Val(:,1),'linewidth',plotWidth);
    end
    if (plot2Val(1,2)>=1e-5)
         plot(plot2Val(:,2),'linewidth',plotWidth);
    end
    if (plot2Val(1,3)>=1e-5)
         plot(plot2Val(:,3),'linewidth',plotWidth);
    end
    plot(ones(1,npts)*ratedVal,'linewidth',plotWidth,'color',ratedCol);
    grid on;
    
    if (commonVal(1,1)>=1e-5)
         plot(commonVal(:,1),'linewidth',plotWidth,'color',[0, 0.5, 0],...
             'linestyle','-.');
    end
    if (commonVal(1,2)>=1e-5)
         plot(commonVal(:,2),'linewidth',plotWidth,'color',[0.75, 0, 0.75],...
         'linestyle','-.');
    end
    if (commonVal(1,3)>=1e-5)
         plot(commonVal(:,3),'linewidth',plotWidth,'color',[0.6350, 0.0780, 0.1840],...
             'linestyle','-.');
    end
    hold off;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am','12:00 pm','12:00 am'});
    
    ylim([floor(min((commonVal), [], 'all'))*0.95 ceil(max([max(plot2Val,[] ,'all'); ratedVal]))*1.0375]) % Max rating are from max current,
                                                                                                                        % normamps and emergamps
    
    title(tileTitle2, 'FontSize', titleFontSize);
    xlim([1 24*nptMult]);
    ax = gca; ax.FontSize = axisFontSize;
    
    xlabel('Time','FontSize', axisFontSize); ylabel('Current (A)','FontSize', axisFontSize);
%     lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southeast','orientation','vertical');
%     lgnd.FontSize = lgndFontSize - decreaseLgnd; clearvars str;

nexttile
hold on;
    if (plot3Val(1,1)>=1e-5)
         plot(plot3Val(:,1),'linewidth',plotWidth);
    end
    if (plot3Val(1,2)>=1e-5)
         plot(plot3Val(:,2),'linewidth',plotWidth);
    end
    if (plot3Val(1,3)>=1e-5)
         plot(plot3Val(:,3),'linewidth',plotWidth);
    end
    plot(ones(1,npts)*ratedVal,'linewidth',plotWidth,'color',ratedCol);
    grid on;
    
    if (commonVal(1,1)>=1e-5)
         plot(commonVal(:,1),'linewidth',plotWidth,'color',[0, 0.5, 0],...
             'linestyle','-.');
    end
    if (commonVal(1,2)>=1e-5)
         plot(commonVal(:,2),'linewidth',plotWidth,'color',[0.75, 0, 0.75],...
         'linestyle','-.');
    end
    if (commonVal(1,3)>=1e-5)
         plot(commonVal(:,3),'linewidth',plotWidth,'color',[0.6350, 0.0780, 0.1840],...
             'linestyle','-.');
         
    end
    hold off;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am','12:00 pm','12:00 am'});
    
    ylim([floor(min((commonVal), [], 'all'))*0.95 ceil(max([max(plot3Val,[] ,'all'); ratedVal]))*1.0375]) % Max rating are from max current,
                                                                                                                        % normamps and emergamps
    lgnd = legend(string(legndArr),'location','southoutside',...
    'Orientation','horizontal',...
    'FontSize',lgndFontSize);
                                                                                                                        
    title(tileTitle3, 'fontsize', titleFontSize);
    xlim([1 24*nptMult]);
    ax = gca; ax.FontSize = axisFontSize;
    
    xlabel('Time','FontSize', axisFontSize); ylabel('Current (A)','FontSize', axisFontSize);
%     lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southeast','orientation','vertical');
%     lgnd.FontSize = lgndFontSize - decreaseLgnd; clearvars str;

nexttile
hold on;
    if (plot4Val(1,1)>=1e-5)
         plot(plot4Val(:,1),'linewidth',plotWidth);
    end
    if (plot4Val(1,2)>=1e-5)
         plot(plot4Val(:,2),'linewidth',plotWidth);
    end
    if (plot4Val(1,3)>=1e-5)
         plot(plot4Val(:,3),'linewidth',plotWidth);
    end
    plot(ones(1,npts)*ratedVal,'linewidth',plotWidth,'color',ratedCol);
   grid on;
   
    if (commonVal(1,1)>=1e-5)
         plot(commonVal(:,1),'linewidth',plotWidth,'color',[0, 0.5, 0],...
             'linestyle','-.');
    end
    if (commonVal(1,2)>=1e-5)
         plot(commonVal(:,2),'linewidth',plotWidth,'color',[0.75, 0, 0.75],...
         'linestyle','-.');
    end
    if (commonVal(1,3)>=1e-5)
         plot(commonVal(:,3),'linewidth',plotWidth,'color',[0.6350, 0.0780, 0.1840],...
             'linestyle','-.');
         
    end
    hold off;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am','12:00 pm','12:00 am'});
    
    ylim([floor(min((commonVal), [], 'all'))*0.95 ceil(max([max(plot4Val,[] ,'all'); ratedVal]))*1.0375]) % Max rating are from max current,
                                                                                                                        % normamps and emergamps
    
    title(tileTitle4, 'fontsize', titleFontSize);
    xlim([1 24*nptMult]);
    ax = gca; ax.FontSize = axisFontSize;
    
    xlabel('Time','FontSize', axisFontSize); ylabel('Current (A)','FontSize', axisFontSize);
%     lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southeast','orientation','vertical');
%     lgnd.FontSize = lgndFontSize - decreaseLgnd; clearvars str;

%% Plotting 634 currents

busNameSet = 4;% busVSet(2);

axisFontSize = 22;
lgndFontSize = 18;
plotWidth = 5;
titleFontSize = 2;
sgTitleFontSize = 26;

figure('Renderer', 'painters', 'Position', [200 -100 1300 900])
    t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
    sgtitle("Current variation for distribution transformer at bus " + busNames(busNameSet)+ newline + " (base kV = " +  num2str(basekVs(busNameSet)) + " kV) in uncoordinated vs coordinated schemes",'fontweight','bold',...
        'fontsize',sgTitleFontSize);
    
nexttile
hold on;
    if (ev40Uncoord.xfm634I(1,1))
         plot(ev40Uncoord.xfm634I(:,1),'linewidth',plotWidth);
    end
    if (ev40Uncoord.xfm634I(1,2))
         plot(ev40Uncoord.xfm634I(:,2),'linewidth',plotWidth);
    end
    if (ev40Uncoord.xfm634I(1,3))
         plot(ev40Uncoord.xfm634I(:,3),'linewidth',plotWidth);
    end
    plot(ones(1,npts)*(1.1*1000/sqrt(3)/0.4),'linewidth',plotWidth);
    grid on;
    %  
    if (normalCond.xfm634I(1,1))
         plot(normalCond.xfm634I(:,1),'linewidth',plotWidth,'color',[0, 0.5, 0]);
    end
    if (normalCond.xfm634I(1,2))
         plot(normalCond.xfm634I(:,2),'linewidth',plotWidth,'color',[0.75, 0, 0.75]);
    end
    if (normalCond.xfm634I(1,3))
         plot(normalCond.xfm634I(:,3),'linewidth',plotWidth,'color',[0.6350, 0.0780, 0.1840]);
    end
    hold off;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am','12:00 pm','12:00 am'});
    
    ylim([floor(min((normalCond.xfm634I), [], 'all'))*0.95 ceil(max([max(ev40Uncoord.xfm634I,[] ,'all'); 1.1*1000/sqrt(3)/0.4]))*1.0375]) % Max rating are from max current,
                                                                                                                        % normamps and emergamps
    
    title(["40% EV penetration" + newline + "with uncoordinated charging"],'fontsize', titleFontSize);
    xlim([1 24*nptMult]);    
    xlabel('Time','fontsize', axisFontSize); ylabel('Current (A)','fontsize', titleFontSize);
    ax = gca; ax.FontSize = axisFontSize;

nexttile
hold on;
    if (ev40Coord.xfm634I(1,1))
         plot(ev40Coord.xfm634I(:,1),'linewidth',plotWidth);
    end
    if (ev40Coord.xfm634I(1,2))
         plot(ev40Coord.xfm634I(:,2),'linewidth',plotWidth);
    end
    if (ev40Coord.xfm634I(1,3))
         plot(ev40Coord.xfm634I(:,3),'linewidth',plotWidth);
    end
    plot(ones(1,npts)*(1.1*1000/sqrt(3)/0.4),'linewidth',plotWidth);
    grid on;
    if (normalCond.xfm634I(1,1))
         plot(normalCond.xfm634I(:,1),'linewidth',plotWidth,'color',[0, 0.5, 0]);
    end
    if (normalCond.xfm634I(1,2))
         plot(normalCond.xfm634I(:,2),'linewidth',plotWidth,'color',[0.75, 0, 0.75]);
    end
    if (normalCond.xfm634I(1,3))
         plot(normalCond.xfm634I(:,3),'linewidth',plotWidth,'color',[0.6350, 0.0780, 0.1840]);
    end
    hold off;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am','12:00 pm','12:00 am'});
    
    ylim([floor(min((normalCond.xfm634I), [], 'all'))*0.95 ceil(max([max(ev40Coord.xfm634I,[] ,'all'); 1.1*1000/sqrt(3)/0.4]))*1.0375]) % Max rating are from max current,
                                                                                                                        % normamps and emergamps
    
    title( "40% EV penetration" + newline + "with coordinated charging", 'fontsize', titleFontSize);
    xlim([1 24*nptMult]);
    
    str = "Normal Ampere" + newline + "rating";
    
    xlabel('Time','fontsize', titleFontSize); ylabel('Current (A)','fontsize', titleFontSize);
    ax = gca; ax.FontSize = axisFontSize;
%     lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southeast','orientation','vertical');
%     lgnd.FontSize = lgndFontSize - decreaseLgnd; clearvars str;

nexttile
hold on;
    if (ev90Uncoord.xfm634I(1,1))
         plot(ev90Uncoord.xfm634I(:,1),'linewidth',plotWidth);
    end
    if (ev90Uncoord.xfm634I(1,2))
         plot(ev90Uncoord.xfm634I(:,2),'linewidth',plotWidth);
    end
    if (ev90Uncoord.xfm634I(1,3))
         plot(ev90Uncoord.xfm634I(:,3),'linewidth',plotWidth);
    end
    plot(ones(1,npts)*(1.1*1000/sqrt(3)/0.4),'linewidth',plotWidth);
    grid on;
    if (normalCond.xfm634I(1,1))
         plot(normalCond.xfm634I(:,1),'linewidth',plotWidth,'color',[0, 0.5, 0]);
    end
    if (normalCond.xfm634I(1,2))
         plot(normalCond.xfm634I(:,2),'linewidth',plotWidth,'color',[0.75, 0, 0.75]);
    end
    if (normalCond.xfm634I(1,3))
         plot(normalCond.xfm634I(:,3),'linewidth',plotWidth,'color',[0.6350, 0.0780, 0.1840]);
    end
    hold off;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am','12:00 pm','12:00 am'});
    
    ylim([floor(min((normalCond.xfm634I), [], 'all'))*0.95 ceil(max([max(ev90Uncoord.xfm634I,[] ,'all'); 1.1*1000/sqrt(3)/0.4]))*1.0375]) % Max rating are from max current,
                                                                                                                        % normamps and emergamps
    
    title( "90% EV penetration" + newline + "with uncoordinated charging", 'fontsize', titleFontSize);
    xlim([1 24*nptMult]);
    
    str = "Normal Ampere" + newline + "rating";
    
    str = "Normal Ampere" + newline + "rating";
%     lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southoutside','orientation','horizontal');
%     lgnd.FontSize = lgndFontSize; clearvars str;
    
    xlabel('Time','fontsize', titleFontSize); ylabel('Current (A)','fontsize', titleFontSize);
    ax = gca; ax.FontSize = axisFontSize;
%     lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southeast','orientation','vertical');
%     lgnd.FontSize = lgndFontSize - decreaseLgnd; clearvars str;

nexttile
hold on;
    if (ev90Coord.xfm634I(1,1))
         plot(ev90Coord.xfm634I(:,1),'linewidth',plotWidth);
    end
    if (ev90Coord.xfm634I(1,2))
         plot(ev90Coord.xfm634I(:,2),'linewidth',plotWidth);
    end
    if (ev90Coord.xfm634I(1,3))
         plot(ev90Coord.xfm634I(:,3),'linewidth',plotWidth);
    end
   plot(ones(1,npts)*(1.1*1e3/sqrt(3)/0.4),'linewidth',plotWidth);
    grid on;
    if (normalCond.xfm634I(1,1))
         plot(normalCond.xfm634I(:,1),'linewidth',plotWidth,'color',[0, 0.5, 0]);
    end
    if (normalCond.xfm634I(1,2))
         plot(normalCond.xfm634I(:,2),'linewidth',plotWidth,'color',[0.75, 0, 0.75]);
    end
    if (normalCond.xfm634I(1,3))
         plot(normalCond.xfm634I(:,3),'linewidth',plotWidth,'color',[0.6350, 0.0780, 0.1840]);
    end
    hold off;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am','12:00 pm','12:00 am'});
    
    ylim([floor(min((normalCond.xfm634I), [], 'all'))*0.95 ceil(max([max(ev90Coord.xfm634I,[] ,'all'); 1.1*1000/sqrt(3)/0.4]))*1.0375]) % Max rating are from max current,
                                                                                                                        % normamps and emergamps
    title( "90% EV penetration" + newline + "with coordinated charging", 'fontsize', titleFontSize);
    xlim([1 24*nptMult]);
    
    str = "Normal Ampere" + newline + "rating";
    
    xlabel('Time','fontsize', titleFontSize); ylabel('Current (A)','fontsize', titleFontSize);
    ax = gca; ax.FontSize = axisFontSize;

%% Plotting 675 currents

busNameSet = busVSet(2);

axisFontSize = 22;
lgndFontSize = 18;
plotWidth = 5;
titleFontSize = 12;
sgTitleFontSize = 22;

figure('Renderer', 'painters', 'Position', [200 -100 1300 900])
    t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
    sgtitle("Current variation for distribution transformer at bus 675 (base kV = " +  num2str(basekVs(busNameSet)) + " kV)",'fontweight','bold',...
        'fontsize',sgTitleFontSize);
    
nexttile
hold on;
    if (ev40Uncoord.xfm675I(1,1))
         plot(ev40Uncoord.xfm675I(:,1),'linewidth',plotWidth);
    end
    if (ev40Uncoord.xfm675I(1,2))
         plot(ev40Uncoord.xfm675I(:,2),'linewidth',plotWidth);
    end
    if (ev40Uncoord.xfm675I(1,3))
         plot(ev40Uncoord.xfm675I(:,3),'linewidth',plotWidth);
    end
    plot(ones(1,npts)*(1.1*1000/sqrt(3)/0.4),'linewidth',plotWidth);
    hold off; grid on;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am','12:00 pm','12:00 am'});
    
    ylim([floor(min((ev40Uncoord.xfm675I), [], 'all'))*0.95 ceil(max([max(ev40Uncoord.xfm675I,[] ,'all'); 1.1*1000/sqrt(3)/0.4]))*1.0375]) % Max rating are from max current,
                                                                                                                        % normamps and emergamps
    
    title('Normal conditions','fontsize', titleFontSize);
    xlim([1 24*nptMult]);    
    xlabel('Time','fontsize', axisFontSize); ylabel('Current (A)','fontsize', titleFontSize);
    ax = gca; ax.FontSize = axisFontSize;

nexttile
hold on;
    if (ev90Uncoord.xfm675I(1,1))
         plot(ev90Uncoord.xfm675I(:,1),'linewidth',plotWidth);
    end
    if (ev90Uncoord.xfm675I(1,2))
         plot(ev90Uncoord.xfm675I(:,2),'linewidth',plotWidth);
    end
    if (ev90Uncoord.xfm675I(1,3))
         plot(ev90Uncoord.xfm675I(:,3),'linewidth',plotWidth);
    end
    plot(ones(1,npts)*(1.1*1000/sqrt(3)/0.4),'linewidth',plotWidth);
    hold off; grid on;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am','12:00 pm','12:00 am'});
    
    ylim([floor(min((ev90Uncoord.xfm675I), [], 'all'))*0.95 ceil(max([max(ev90Uncoord.xfm675I,[] ,'all'); 1.1*1000/sqrt(3)/0.4]))*1.0375]) % Max rating are from max current,
                                                                                                                        % normamps and emergamps
    
    title('20 % EV penetration', 'fontsize', titleFontSize);
    xlim([1 24*nptMult]);
    
    str = "Normal Ampere" + newline + "rating";
    
    xlabel('Time','fontsize', titleFontSize); ylabel('Current (A)','fontsize', titleFontSize);
    ax = gca; ax.FontSize = axisFontSize;
%     lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southeast','orientation','vertical');
%     lgnd.FontSize = lgndFontSize - decreaseLgnd; clearvars str;

nexttile
hold on;
    if (ev90Coord.xfm675I(1,1))
         plot(ev90Coord.xfm675I(:,1),'linewidth',plotWidth);
    end
    if (ev90Coord.xfm675I(1,2))
         plot(ev90Coord.xfm675I(:,2),'linewidth',plotWidth);
    end
    if (ev90Coord.xfm675I(1,3))
         plot(ev90Coord.xfm675I(:,3),'linewidth',plotWidth);
    end
    plot(ones(1,npts)*(1.1*1000/sqrt(3)/0.4),'linewidth',plotWidth);
    hold off; grid on;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am','12:00 pm','12:00 am'});
    
    ylim([floor(min((ev90Coord.xfm675I), [], 'all'))*0.95 ceil(max([max(ev90Coord.xfm675I,[] ,'all'); 1.1*1000/sqrt(3)/0.4]))*1.0375]) % Max rating are from max current,
                                                                                                                        % normamps and emergamps
    
    title('40 % EV penetration', 'fontsize', titleFontSize);
    xlim([1 24*nptMult]);
    
    str = "Normal Ampere" + newline + "rating";
    
    str = "Normal Ampere" + newline + "rating";
    lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southoutside','orientation','horizontal');
    lgnd.FontSize = lgndFontSize; clearvars str;
    
    xlabel('Time','fontsize', titleFontSize); ylabel('Current (A)','fontsize', titleFontSize);
    ax = gca; ax.FontSize = axisFontSize;
%     lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southeast','orientation','vertical');
%     lgnd.FontSize = lgndFontSize - decreaseLgnd; clearvars str;

nexttile
hold on;
    if (ev90Coord.xfm675I(1,1))
         plot(ev90Coord.xfm675I(:,1),'linewidth',plotWidth);
    end
    if (ev90Coord.xfm675I(1,2))
         plot(ev90Coord.xfm675I(:,2),'linewidth',plotWidth);
    end
    if (ev90Coord.xfm675I(1,3))
         plot(ev90Coord.xfm675I(:,3),'linewidth',plotWidth);
    end
   plot(ones(1,npts)*(1.1*1e3/sqrt(3)/0.4),'linewidth',plotWidth);
    hold off; grid on;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am','12:00 pm','12:00 am'});
    
    ylim([floor(min((ev90Coord.xfm675I), [], 'all'))*0.95 ceil(max([max(ev90Coord.xfm675I,[] ,'all'); 1.1*1000/sqrt(3)/0.4]))*1.0375]) % Max rating are from max current,
                                                                                                                        % normamps and emergamps
    
    title('90 % EV penetration', 'fontsize', titleFontSize);
    xlim([1 24*nptMult]);
    
    str = "Normal Ampere" + newline + "rating";
    
    xlabel('Time','fontsize', titleFontSize); ylabel('Current (A)','fontsize', titleFontSize);
    ax = gca; ax.FontSize = axisFontSize;
%     lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southeast','orientation','vertical');
%     lgnd.FontSize = lgndFontSize - decreaseLgnd; clearvars str;
    
% nexttile
% hold on;
%     if (normal.xfm675I(1,1))
%          plot(normal.xfm675I(:,1),'linewidth',plotWidth);
%     end
%     if (normal.xfm675I(1,2))
%          plot(normal.xfm675I(:,2),'linewidth',plotWidth);
%     end
%     if (normal.xfm675I(1,3))
%          plot(normal.xfm675I(:,3),'linewidth',plotWidth);
%     end
%     plot(ones(1,npts)*(1.1*1e3/sqrt(3)/0.4),'linewidth',plotWidth);
%     hold off; grid on;
%     
%     xticks(linspace(1,24*nptMult,5));
%     xticklabels({'12:00 am','6:00 am'...
%         '12:00 pm','6:00 pm','12:00 am'});
%     
%     ylim([floor(min((normal.xfm675I), [], 'all'))/500 ceil(max([max(normal.xfm675I,[] ,'all'); 1.1*1e3/sqrt(3)/0.4]))*1.075]) % Max rating are from max current,
%                                                                                                                         % normamps and emergamps
%     
%     title('Normal conditions');
%     xlim([1 24*nptMult]);
%     ax = gca; ax.FontSize = axisFontSize - decreaseAxis;
%     
%     str = "Normal Ampere" + newline + "rating";
%     lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southeast','orientation','vertical');
%     lgnd.FontSize = lgndFontSize - decreaseLgnd; clearvars str;
% 
% nexttile
% hold on;
%     if (ev20.xfm675I(1,1))
%          plot(ev20.xfm675I(:,1),'linewidth',plotWidth);
%     end
%     if (ev20.xfm675I(1,2))
%          plot(ev20.xfm675I(:,2),'linewidth',plotWidth);
%     end
%     if (ev20.xfm675I(1,3))
%          plot(ev20.xfm675I(:,3),'linewidth',plotWidth);
%     end
%     plot(ones(1,npts)*(1.1*1e3/sqrt(3)/0.4),'linewidth',plotWidth);
%     hold off; grid on;
%     
%     xticks(linspace(1,24*nptMult,5));
%     xticklabels({'12:00 am','6:00 am'...
%         '12:00 pm','6:00 pm','12:00 am'});
%     
%     ylim([floor(min((ev20.xfm675I), [], 'all'))/500 ceil(max([max(ev20.xfm675I,[] ,'all'); 1.1*1e3/sqrt(3)/0.4]))*1.075]) % Max rating are from max current,
%                                                                                                                         % normamps and emergamps                                                                                                                        % normamps and emergamps
%     
%     title('20 % EV penetration', 'fontsize', titleFontSize - decreaseTitle);
%     xlim([1 24*nptMult]);
%     ax = gca; ax.FontSize = axisFontSize - decreaseAxis;
%     
%     str = "Normal Ampere" + newline + "rating";
% %     lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southeast','orientation','vertical');
% %     lgnd.FontSize = lgndFontSize - decreaseLgnd; clearvars str;
% 
% nexttile
% hold on;
%     if (ev40.xfm675I(1,1))
%          plot(ev40.xfm675I(:,1),'linewidth',plotWidth);
%     end
%     if (ev40.xfm675I(1,2))
%          plot(ev40.xfm675I(:,2),'linewidth',plotWidth);
%     end
%     if (ev40.xfm675I(1,3))
%          plot(ev40.xfm675I(:,3),'linewidth',plotWidth);
%     end
%     plot(ones(1,npts)*(1.1*1e3/sqrt(3)/0.4),'linewidth',plotWidth);
%     hold off; grid on;
%     
%     xticks(linspace(1,24*nptMult,5));
%     xticklabels({'12:00 am','6:00 am'...
%         '12:00 pm','6:00 pm','12:00 am'});
%     
%     ylim([floor(min((ev40.xfm675I), [], 'all'))/500 ceil(max([max(ev40.xfm675I,[] ,'all'); 1.1*1e3/sqrt(3)/0.4]))*1.075]) % Max rating are from max current,
%                                                                                                                         % normamps and emergamps                                                                                                                        % normamps and emergamps
%     
%     title('40 % EV penetration', 'fontsize', titleFontSize - decreaseTitle);
%     xlim([1 24*nptMult]);
%     ax = gca; ax.FontSize = axisFontSize - decreaseAxis;
%     
%     str = "Normal Ampere" + newline + "rating";
% %     lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southeast','orientation','vertical');
% %     lgnd.FontSize = lgndFontSize - decreaseLgnd; clearvars str;
% 
% nexttile
% hold on;
%     if (ev90.xfm675I(1,1))
%          plot(ev90.xfm675I(:,1),'linewidth',plotWidth);
%     end
%     if (ev90.xfm675I(1,2))
%          plot(ev90.xfm675I(:,2),'linewidth',plotWidth);
%     end
%     if (ev90.xfm675I(1,3))
%          plot(ev90.xfm675I(:,3),'linewidth',plotWidth);
%     end
%     plot(ones(1,npts)*(1.1*1e3/sqrt(3)/0.4),'linewidth',plotWidth);
%     hold off; grid on;
%     
%     xticks(linspace(1,24*nptMult,5));
%     xticklabels({'12:00 am','6:00 am'...
%         '12:00 pm','6:00 pm','12:00 am'});
%     
%     ylim([floor(min((ev90.xfm675I), [], 'all'))/500 ceil(max([max(ev90.xfm675I,[] ,'all'); 1.1*1e3/sqrt(3)/0.4]))*1.075]) % Max rating are from max current,
%                                                                                                                         % normamps and emergamps                                                                                                                        % normamps and emergamps
%     
%     title('90 % EV penetration', 'fontsize', titleFontSize - decreaseTitle);
%     xlim([1 24*nptMult]);
%     ax = gca; ax.FontSize = axisFontSize - decreaseAxis;
%     
%     str = "Normal Ampere" + newline + "rating";
%     lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southeast','orientation','vertical');
%     lgnd.FontSize = lgndFontSize - decreaseLgnd; clearvars str;

%% Plotting 652 currents

busNameSet = busVSet(3);

axisFontSize = 22;
lgndFontSize = 18;
plotWidth = 5;
titleFontSize = 12;
sgTitleFontSize = 22;

figure('Renderer', 'painters', 'Position', [200 -100 1300 900])
    t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
    sgtitle("Current variation for distribution transformer at bus 652 (base kV = " +  num2str(basekVs(busNameSet)) + " kV)",'fontweight','bold',...
        'fontsize',sgTitleFontSize);
    
nexttile
hold on;
    if (ev40Uncoord.xfm652I(1,1)>=1e-5)
         plot(ev40Uncoord.xfm652I(:,1),'linewidth',plotWidth);
    end
    if (ev40Uncoord.xfm652I(1,2)>=1e-5)
         plot(ev40Uncoord.xfm652I(:,2),'linewidth',plotWidth);
    end
    if (ev40Uncoord.xfm652I(1,3)>=1e-5)
         plot(ev40Uncoord.xfm652I(:,3),'linewidth',plotWidth);
    end
    plot(ones(1,npts)*(1.1*1000/sqrt(3)/0.4),'linewidth',plotWidth,'color',[0.4940, 0.1840, 0.5560]);
    hold off; grid on;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am',...
        '12:00 pm','12:00 am'});
    title('Normal conditions');
    
    ylim([floor(min((ev40Uncoord.xfm652I(:,1)), [], 'all'))*0.95 ceil(max([max(ev40Uncoord.xfm652I,[] ,'all'); 1.1*1000/sqrt(3)/0.4]))*1.0375]) % Max rating are from max current,
                                                                                                                        % normamps and emergamps
    xlabel('Time'); ylabel('Current (A)');
    xlim([1 24*nptMult]);
    ax = gca; ax.FontSize = axisFontSize;

nexttile
hold on;
    if (ev90Uncoord.xfm652I(1,1)>=1e-5)
         plot(ev90Uncoord.xfm652I(:,1),'linewidth',plotWidth);
    end
    if (ev90Uncoord.xfm652I(1,2)>=1e-5)
         plot(ev90Uncoord.xfm652I(:,2),'linewidth',plotWidth);
    end
    if (ev90Uncoord.xfm652I(1,3)>=1e-5)
         plot(ev90Uncoord.xfm652I(:,3),'linewidth',plotWidth);
    end
    plot(ones(1,npts)*(1.1*1000/sqrt(3)/0.4),'linewidth',plotWidth,'color',[0.4940, 0.1840, 0.5560]);
    hold off; grid on;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am',...
        '12:00 pm','12:00 am'});
    
    ylim([floor(min((ev90Uncoord.xfm652I(:,1)), [], 'all'))*0.95 ceil(max([max(ev90Uncoord.xfm652I,[] ,'all'); 1.1*1000/sqrt(3)/0.4]))*1.0375]) % Max rating are from max current,
                                                                                                                        % normamps and emergamps
    
    title('20 % EV penetration', 'fontsize', titleFontSize);
    xlabel('Time'); ylabel('Current (A)');
    xlim([1 24*nptMult]);
    ax = gca; ax.FontSize = axisFontSize;
    
    str = "Normal Ampere" + newline + "rating";
    
%     lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southeast','orientation','vertical');
%     lgnd.FontSize = lgndFontSize - decreaseLgnd; clearvars str;

nexttile
hold on;
    if (ev90UnCoord.xfm652I(1,1)>=1e-5)
         plot(ev90UnCoord.xfm652I(:,1),'linewidth',plotWidth);
    end
    if (ev90UnCoord.xfm652I(1,2)>=1e-5)
         plot(ev90UnCoord.xfm652I(:,2),'linewidth',plotWidth);
    end
    if (ev90UnCoord.xfm652I(1,3)>=1e-5)
         plot(ev90UnCoord.xfm652I(:,3),'linewidth',plotWidth);
    end
    plot(ones(1,npts)*(1.1*1000/sqrt(3)/0.4),'linewidth',plotWidth,'color',[0.4940, 0.1840, 0.5560]);
    hold off; grid on;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am',...
        '12:00 pm','12:00 am'});
    
    ylim([floor(min((ev90UnCoord.xfm652I(:,1)), [], 'all'))*0.95 ceil(max([max(ev90UnCoord.xfm652I,[] ,'all'); 1.1*1000/sqrt(3)/0.4]))*1.0375]) % Max rating are from max current,
                                                                                                                        % normamps and emergamps
    
    title('40 % EV penetration', 'fontsize', titleFontSize);
    xlabel('Time'); ylabel('Current (A)');
    xlim([1 24*nptMult]);
    ax = gca; ax.FontSize = axisFontSize;
    
    str = "Normal Ampere" + newline + "rating";
    
    str = "Normal Ampere" + newline + "rating";
    lgnd = legend({'Phase a',str},'Location','southoutside','orientation','horizontal');
    lgnd.FontSize = lgndFontSize; clearvars str;
    
%     lgnd = legend({'Phase a','Phase b','Phase c',str},'Location','southeast','orientation','vertical');
%     lgnd.FontSize = lgndFontSize - decreaseLgnd; clearvars str;

nexttile
hold on;
    if (ev90Coord.xfm652I(1,1)>=1e-5)
         plot(ev90Coord.xfm652I(:,1),'linewidth',plotWidth);
    end
    if (ev90Coord.xfm652I(1,2)>=1e-5)
         plot(ev90Coord.xfm652I(:,2),'linewidth',plotWidth);
    end
    if (ev90Coord.xfm652I(1,3)>=1e-5)
         plot(ev90Coord.xfm652I(:,3),'linewidth',plotWidth);
    end
   plot(ones(1,npts)*(1.1*1e3/sqrt(3)/0.4),'linewidth',plotWidth,'color',[0.4940, 0.1840, 0.5560]);
    hold off; grid on;
    
    xticks(linspace(1,24*nptMult,3));
    xticklabels({'12:00 am',...
        '12:00 pm','12:00 am'});
    
    ylim([floor(min((ev90Coord.xfm652I(:,1)), [], 'all'))*0.95 ceil(max([max(ev90Coord.xfm652I,[] ,'all'); 1.1*1000/sqrt(3)/0.4]))*1.0375]) % Max rating are from max current,
                                                                                                                        % normamps and emergamps
    
    title('90 % EV penetration', 'fontsize', titleFontSize);
    xlabel('Time'); ylabel('Current (A)');
    xlim([1 24*nptMult]);
    ax = gca; ax.FontSize = axisFontSize;
    
    str = "Normal Ampere" + newline + "rating";

%% Plotting 24-hr power flow
close all;

axisFontSize = 22;
lgndFontSize = 10;
plotWidth = 5;
titleFontSize = 12;

plot1Val = ev40Uncoord.xfmkw675;
plot2Val = ev90Uncoord.xfmkw675;
plot3Val = ev40Coord.xfmkw675;
plot4Val = ev90Coord.xfmkw675;
commonVal = normalCond.xfmkw675;

ratedVal = ratedDisTr675kw;

titleText = "Bus 675 distribution transformer active power flow" + newline + ...
    "in uncoordinated vs coordinated schemes";

figure('Renderer', 'painters', 'Position', [200 -100 1100 700])
t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile;
hold on;
plot(sum(commonVal, 2),'linewidth',plotWidth,'color',[0    0.6510    0.1569]);
plot(sum(plot1Val, 2),'linewidth',plotWidth,'color','red','linestyle','-.');
plot(sum(plot2Val, 2),'linewidth',plotWidth,'color','blue','linestyle','-.');
plot(sum(plot3Val, 2),'linewidth',plotWidth,'color','red');
plot(sum(plot4Val, 2),'linewidth',plotWidth,'color','blue');
plot(ones(1,npts)*ratedVal,'linewidth',plotWidth,'color',[125, 50, 168]/255);
hold off;

xticks(linspace(1,24*nptMult,7));
    xticklabels({'12:00 am','04:00 am','08:00 am',...
        '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
grid on;
title(titleText, 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
str0 = "Normal conditions";
str1 = "40% EV penetration" + newline + "with uncoordinated charging";
str2 = "90% EV penetration" + newline + "with uncoordinated charging";
str3 = "40% EV penetration" + newline + "with coordinated charging";
str4 = "90% EV penetration" + newline + "with coordinated charging";
str5 = "Maximum nominal demand" + newline + "under normal conditions";
lgnd = legend({str0,str1,str2,str3,str4,str5},...
    'Location','northwest','orientation','vertical');
lgnd.FontSize = lgndFontSize;
ylim([floor(min(sum(commonVal, 2)))*0.95 ((max(sum(plot2Val, 2))/500))*500*1.05])
xlim([1 24*nptMult]);

%% Losses as a percentage

lossesA = [(trapz(sum(ev40Uncoord.subkw, 2)) - trapz(sum(ev40Uncoord.resLdkw, 2)...
    + sum(ev40Uncoord.industLdkw, 2)))/trapz(sum(ev40Uncoord.subkw, 2))*100;
    (trapz(sum(ev90Uncoord.subkw, 2)) - trapz(sum(ev90Uncoord.resLdkw, 2)...
    + sum(ev90Uncoord.industLdkw, 2)) - trapz(sum(ev90Uncoord.PevTot, 2)))...
    /trapz(sum(ev40Uncoord.subkw, 2))*100;
    (trapz(sum(ev40Coord.subkw, 2)) - trapz(sum(ev40Coord.resLdkw, 2)...
    + sum(ev40Coord.industLdkw, 2)) - trapz(sum(ev40Coord.PevTot, 2)))...
    /trapz(sum(ev40Uncoord.subkw, 2))*100;
    (trapz(sum(ev90Coord.subkw, 2)) - trapz(sum(ev90Coord.resLdkw, 2)...
    + sum(ev90Coord.industLdkw, 2)) - trapz(sum(ev90Coord.PevTot, 2)))...
    /trapz(sum(ev40Uncoord.subkw, 2))*100];% , [1:npts]/nptMult);

losses = [sum(normalCond.PLoss, 1);sum(ev40Uncoord.PLoss, 1); sum(ev90Uncoord.PLoss, 1); ...
    sum(ev40Coord.PLoss, 1); sum(ev90Coord.PLoss, 1)];

% different scenario losses compared with normal input
temp1 = sum(losses(1,:), 'all')./sum(ev40Uncoord.subkw,'all')*100;
temp2 = sum(losses(2,:), 'all')./sum(ev40Uncoord.subkw,'all')*100;
temp3 = sum(losses(3,:), 'all')./sum(ev40Uncoord.subkw,'all')*100;
temp4 = sum(losses(4,:), 'all')./sum(ev40Uncoord.subkw,'all')*100;

% % different scenario losses compared with their respective input
% Loss2 = [sum(losses(1,:), 'all')./sum(normal.subkw,'all')*100;
% sum(losses(2,:), 'all')./sum(ev20.subkw,'all')*100;
% sum(losses(3,:), 'all')./sum(ev40.subkw,'all')*100;
% sum(losses(4,:), 'all')./sum(ev90.subkw,'all')*100];

Loss2 = [sum(losses(1,:), 'all');
sum(losses(2,:), 'all');
sum(losses(3,:), 'all');
sum(losses(4,:), 'all')];

Loss = [temp1 temp2 temp3 temp4];
increaseLoss = [Loss(1) (Loss(2:4)-Loss(1))*100/Loss(1)];

increaseLoss2 = [Loss2(1); (Loss2(2:4)-Loss2(1))*100/Loss2(1)];

axisFontSize = 22;
lgndFontSize = 16;
plotWidth = 5;
titleFontSize = 12;

figure('Renderer', 'painters', 'Position', [200 -100 1100 700])
t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile;
plot(losses(1,:)','linewidth',plotWidth,'color',[0    0.6510    0.1569]); hold on;
plot(losses(2,:)','linewidth',plotWidth-0.4,'color','red','linestyle','-.');
plot(losses(3,:)','linewidth',plotWidth-0.8, 'color','blue','linestyle','-.');
plot(losses(4,:)','linewidth',plotWidth-1.2, 'color', 'red');
plot(losses(5,:)','linewidth',plotWidth-1.6, 'color','blue'); hold off;
ylim([floor(min(min(losses,[],'all'))*0.8) ceil(max(losses,[],'all'))+5])
xlim([1 24*nptMult]); grid on;
xticks(linspace(1,24*nptMult,7));
xticklabels({'12:00 am','04:00 am','08:00 am',...
    '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
xlabel('Time'); ylabel('Varying losses (kW)');
title("Total daily losses of the system in uncoordinated vs" + newline + "coordinated charging", 'FontSize', titleFontSize); xlabel('Time',...
    'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
ax = gca; ax.FontSize = axisFontSize;
% str = "Sum of" + newline + "all phases";
lgnd = legend({"Normal conditions","40% EV penetration with" + newline + "uncoordinated charging","90% EV penetration with" + newline + "uncoordinated charging",...
    "40% EV penetration with" + newline + "coordinated charging","90% EV penetration with" + newline + "coordinated charging"},'Location',...
    'northwest','orientation','vertical');
lgnd.FontSize = lgndFontSize; clearvars str;