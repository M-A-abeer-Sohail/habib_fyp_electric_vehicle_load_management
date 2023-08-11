%%----------------- main.m junk code ------------------------
% %----------------Junk Code---------------------
% 
% % xfmkw634(xfmkw634==ISNAN) = 0;
% % xfmkvar634(xfmkvar634==ISNAN) = 0;
% % 
% % xfmkw675(xfmkw675==ISNAN) = 0;
% % xfmkvar675(xfmkvar675==ISNAN) = 0;
% % 
% % xfmkw611(xfmkw611==ISNAN) = 0;
% % xfmkvar611(xfmkvar611==ISNAN) = 0;
% % 
% % xfmkw652(xfmkw652==ISNAN) = 0;
% % xfmkvar652(xfmkvar652==ISNAN) = 0;
% % 
% % xfmkw646(xfmkw646==ISNAN) = 0;
% % xfmkvar646(xfmkvar646==ISNAN) = 0;
% 
% % plot(5500*1.1/11/sqrt(3)*ones(1,96))
% % plot(5500*1.5/11/sqrt(3)*ones(1,96))
% % figure; plot(subI)
% % hold on
% % plot(5500*1.1/11/sqrt(3)*ones(1,96))
% % plot(5500*1.5/11/sqrt(3)*ones(1,96))
% 
% %%
% % This will contain all the losses
% lossesAll = [lossesAll (sum(subkw, 2) - sum(resLdkw, 2) - sum(industLdkw, 2))];
% 
% % lossesPercent = [lossesPercent sum(lossesAll(:,1))*100/sum(sum(subkw, 2))];
% 
% %% 20% EV penetration set up
% 
% close all;
% 
% %% EV car numbers
% 
% pkCar = [224000 233000 244000];
% khiCar = 0.2*pkCar;
% evCarProp = [0.2 0.4 0.9];
% evCarKhi = evCarProp.*khiCar;
% % evCar13 = [ceil(0.01*evCarKhi(1)) floor(0.01*evCarKhi(2:3)); % This is the number of EV cars in the IEEE 13 bus
% evCar13 = [90 186 439];
% 
% pkBike = [2300000 2800000 3900000];
% khiBike = 0.2*pkBike;
% evBikeProp = [0.1 0.5 0.9];
% evBikeKhi = evBikeProp.*khiBike;
% % evBike13 = floor(0.01*evBikeKhi); % This is the number of EV cars in the IEEE 13 bus
% evBike13 = [460 2800 7000];
% 
% %%
% numEVLds = 13;
% 
% evBusNames = [634, 671, 675, 692, 680]; % In order of increasing load
% 
% Vi=ones(1,numEVLds); % Initial bus reference voltage, but I think I should
%                      % lay a condition that relates to Haidar paper
%                      % These are voltages per phase
% 
% % Will need to fix these...
% evAtBus = readmatrix('data/evVeh20.csv');
% carAtBus=evAtBus(:,1:end-2);
% bikeAtBus=evAtBus(:,end-1:end);
% 
% evBikeV = [36 48]; % in volts
% evBikeAh = [20 20]; % in Ah
% evBikeChTime = [7 7]; % in hrs
% 
% % except EV volts, data taken from https://www.mobilityhouse.com/int_en/knowledge-center/charging-time-summary
% % check battery voltages.xls
% evCarV = [356 240 360 375 356 240 346 330 300 360 350]; % in Volts
% evCarkwh = [64 28 27 30 64 4.4 8.8 16 12 24 30]; % in kwh
% evCarChTime = [4.5 4.5 6 7.5 9.5 1.5 2.5 6 5 5.5 7]; % in hrs, wallbox
% 
% evPowerFactor = 0.9;
% 
% busNumPh = [3 3 3 2 1 1];
% 
% str_part2 = ["New Load.EV_634a Bus1=634.1 Phases=1 Conn=Wye  Model=1 kV=0.23 vminpu=0.3 ",...
%             'New Load.EV_634b Bus1=634.2 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
%             'New Load.EV_634c Bus1=634.3 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
%             'New Load.EV_675a Bus1=675.1 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
%             'New Load.EV_675b Bus1=675.2 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
%             'New Load.EV_675c Bus1=675.3 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
%             'New Load.EV_680a Bus1=680.1 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
%             'New Load.EV_680b Bus1=680.2 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
%             'New Load.EV_680c Bus1=680.3 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
%             "New Load.EV_646a Bus1=646.1 Phases=1 Conn=Wye  Model=1 kV=0.23 vminpu=0.3 ",...
%             "New Load.EV_646b Bus1=646.2 Phases=1 Conn=Wye  Model=1 kV=0.23 vminpu=0.3 ",...
%             "New Load.EV_652a Bus1=652.1 Phases=1 Conn=Wye  Model=1 kV=0.23 vminpu=0.3 ",...
%             "New Load.EV_611c Bus1=611.3 Phases=1 Conn=Wye  Model=1 kV=0.23 vminpu=0.3 "];
% 
% % peakTime = 12;
% carPdfMult = 1;
% 
% carPdf = normpdf(linspace(-1.1,0.5,npts),0.2,0.13) + 0.4*rand(1,npts); % this will multiply
%                                                                        % with the
%                                                                        % EV load to
%                                                                        % simulate gradual rise,
%                                                                        % also with some noise
% carPdf = 0.05+carPdfMult*carPdf/max(carPdf);% normalizing, 
%                                             % and scaling 
%                                             % up/down with the mult
% plot(carPdf);
% 
% %%
% % % Arrays that'll hold per-phase voltage values
% % % These will be cleared for each case
% % V1pu = [];
% V2pu = [];
% V3pu = [];
% 
% % power flows for substation and distribution transformer
% % These will be cleared for each case
% subkw = [];
% subkvar = [];
% xfmkw = [];
% xfmkvar = [];
% PLoss = []; 
% QLoss = [];
% PevAll = [];
% QevAll = [];
% Pev = [];
% Qev = [];
% 
% ViSetting = [4 6 11 10 7];
% % I confirmed the voltages and their buses
% % This is correct as per str_part2
% 
% %----------------Junk Code-------------------
% 
% % Generating random EV's
% % evAtBus = randi([5,75],[numEVLds, 85]);
% 
% % IMPORTANT DATA: This shows number of EV's at a bus
% % evAtBus = [20 10 30 15 15];
% 
% % The following loads the number of EV car models at each bus, as a 5 x 85 matrix
% 
% %% 20% EV penetration all plots
% 
% for i = 1:npts
%     kwLoopRes = reskw*resMult(i); kvarLoopRes = reskvar*resMult(i); % The particular power
%                                                                         % values for the i-th hour
%     kwLoopInd = industkw*industMult(i); kvarLoopInd = industkvar*industMult(i);
%     
%     [Pev, Qev] = ev_load(Vi, carAtBus*carPdf(i));
% 
% %     Pev = rowCopyX(Pev, 3);
% %     Qev = rowCopyX(Qev, 3);
% 
%     makingLoads([str_part, str_part2], [kwLoopInd(1) kwLoopRes(1:6) kwLoopInd(2:end) kwLoopRes(7:end),...
%         Pev], [kvarLoopInd(1) kvarLoopRes(1:6) kvarLoopInd(2:end) kvarLoopRes(7:end),...
%         Qev]); % make the loadss.txt file
%     
%     PevAll = [PevAll; Pev];
%     QevAll = [QevAll; Qev];
%     
%     DailyLoadFlow; % compute the "daily" load flow from the files in OpenDSS
%     V1pu = [V1pu; Vpu1hr(:,1)']; % Get the per-phase voltage values for phase a
%     V2pu = [V2pu; Vpu1hr(:,2)']; % Get the per-phase voltage values for phase b
%     V3pu = [V3pu; Vpu1hr(:,3)']; % Get the per-phase voltage values for phase c
%     tempP = xlsread('MasterIEEE13_EXP_P_BYPHASE.CSV',1,'E2:J3'); % Read the powers per phase from the file
%     [tempsubkw, tempsubkvar, tempxfmkw, tempxfmkvar] = getPs(tempP); % Extract the power
%                                                                      % values for substation
%                                                                      % and distr trafo.
%     subkw = [subkw; tempsubkw]; % Store the kW power flow for substation
%     subkvar = [subkvar; tempsubkvar]; % Store the kVAR power flow for substation
%     xfmkw = [xfmkw; tempxfmkw]; % Store the kW power flow for dist. trafo.
%     xfmkvar = [xfmkvar; tempxfmkvar]; % Store the kVAR power flow for dist. trafo.
%     
% %     for k = 1:numEVBus
% %         Vi(k)=V1pu(i,ViSetting(k));
% %     end
%     
%     % Putting in the phase of the load
%     Vi(1) = V1pu(i,ViSetting(1));
%     Vi(2) = V2pu(i,ViSetting(1));
%     Vi(3) = V3pu(i,ViSetting(1));
%     
%     Vi(4) = V1pu(i,ViSetting(2));
%     Vi(5) = V2pu(i,ViSetting(2));
%     Vi(6) = V3pu(i,ViSetting(2));
%     
%     Vi(7) = V1pu(i,ViSetting(3));
%     Vi(8) = V2pu(i,ViSetting(3));
%     Vi(9) = V3pu(i,ViSetting(3));
%     
%     Vi(10) = V1pu(i,ViSetting(4));
%     Vi(11) = V2pu(i,ViSetting(4));
%     Vi(12) = V3pu(i,ViSetting(4));
%     
%     Vi(13) = V1pu(i,ViSetting(5));
%     Vi(14) = V2pu(i,ViSetting(5));
%     Vi(15) = V3pu(i,ViSetting(5));
% 
%     [PLossTemp, QLossTemp] = getLinePower(readtable("MasterIEEE13_EXP_LOSSES.CSV"));
%     PLoss = [PLoss; PLossTemp];
%     QLoss = [QLoss; QLossTemp];
%     
% end
% 
% V1pu(V1pu==0) = NaN;
% V2pu(V2pu==0) = NaN;
% V3pu(V3pu==0) = NaN;
% 
% save('ev20Pent','V1pu','V2pu','V3pu',...
%     'ratedDistTrafokw','xfmkw', 'xfmkvar',...
%     'subkw','subkvar', 'ratedSubkw', 'houseLdkw','industLdkw',...
%     'houseLdkvar','industLdkvar','lossesBus','efficiencySub','efficiencyDistTrafo',...
%     'PevAll', 'QevAll');
% 
% % Load variation for all loads
% figure;
% 
% plot(sum(resLdkw, 2) + sum(industLdkw, 2), 'linewidth', plotWidth, 'Color','blue','Marker','o','MarkerFaceColor','white');  hold on
% plot(sum(resLdkw, 2) + sum(industLdkw, 2) + sum(PevAll, 2),'Color','red','Marker','o','MarkerFaceColor','white',...
%     'linewidth',plotWidth);
% grid on; xticks(linspace(1,24*nptMult,7));
% xticklabels({'12:00 am','04:00 am','08:00 am','12:00 pm','04:00 pm','08:00 pm','12:00 am'});
% title('Load variation with 20% electric vehicle penetration', 'FontSize', titleFontSize); xlabel('Time',...
%     'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
% ax = gca; ax.FontSize = axisFontSize;
% lgnd = legend(["Normal scenario" "Load with 20% EV penetration"],'Location','southoutside','orientation','horizontal');
% lgnd.FontSize = lgndFontSize;
% xlim([1 24*nptMult]);
% ylim([(floor(min(sum(resLdkw, 2) + sum(industLdkw, 2))/500)*500 - 500),...
%     (ceil(max(sum(resLdkw, 2) + sum(industLdkw, 2)+sum(PevAll, 2))/500)*500 + 500)])
% 
% % Voltages (check busNames for voltage values)
% for i =  2:length(busNames) % ignoring SOURCEBUS
%     figure('Renderer', 'painters', 'Position', [200 100 1100 700]);
%     hold on;
%     if ~isnan(V1pu(1,i))
%          plot(V1pu(:,i),'linewidth',2);
%     end
%     if ~isnan(V2pu(1,i))
%          plot(V2pu(:,i),'linewidth',2);
%     end
%     if ~isnan(V3pu(1,i))
%          plot(V3pu(:,i),'linewidth',2);
%     end
%     plot(1.05*ones(1,npts),'linewidth',2,'linestyle','--','color','black');
%     plot(0.95*ones(1,npts),'linewidth',2,'linestyle','--','color','black');
%     hold off
%     
%     xlabel('Time in hours', 'FontSize', axisFontSize);ylabel('Voltage in PU', 'FontSize', axisFontSize);
%     title(strcat("Voltage variation for bus", " ", busNames(i),", ","with 20% EV penetration (base kV_{\phi} ="," ",...
%         num2str(round(basekVs(i)/sqrt(3),3))," ","kV_{rms})"),'FontSize', titleFontSize);
%     
%     legndArr = {};
%     if ~isnan(V1pu(1,i))
%         legndArr(end+1) = {"Phase a"};
%     end
%     if ~isnan(V2pu(1,i))
%          legndArr(end+1) = {"Phase b"};
%     end
%     if ~isnan(V3pu(1,i))
%          legndArr(end+1) = {"Phase c"};
%     end
%     str1 = "1.05 p.u." + newline + "limit";
%     str2 = "0.95 p.u." + newline + "limit";
%     legndArr(end+1) = {str1};
%     legndArr(end+1) = {str2};
%     clearvars str1 str2
%     lgnd = legend(string(legndArr),'Location','southoutside','orientation','horizontal');
%     lgnd.FontSize = lgndFontSize; clear str;
%     grid on;
%     ax = gca; ax.FontSize = axisFontSize;
%     
%     xticks(linspace(1,24*nptMult,7));
%     xticklabels({'12:00 am','04:00 am','08:00 am','12:00 pm','04:00 pm','08:00 pm','12:00 am'});
%     ylim([min([V1pu(:,i)' V2pu(:,i)' V3pu(:,i)'])-0.01 1.1])
%     xlim([1 24*nptMult]);
% end
% 
% % Plotting 24-hr flow at substation
% figure('Renderer', 'painters', 'Position', [200 100 1100 700]);
% plot(subkw(:,1),'linewidth',2); hold on;
% plot(subkw(:,2),'linewidth',2); 
% plot(subkw(:,3),'linewidth',2);
% plot(sum(subkw, 2),'linewidth',2);
% plot(ones(1,npts)*ratedSubkw,'linewidth',2);
% hold off;
% 
% xticks(linspace(1,24*nptMult,7));
%     xticklabels({'12:00 am','04:00 am','08:00 am',...
%         '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
% grid on;
% title('Substation (bus 650) active power flow with 20% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
%     'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
% ax = gca; ax.FontSize = axisFontSize;
% str = "Sum of" + newline + "all phases";
% lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
% lgnd.FontSize = lgndFontSize; clearvars str;
% grid on;
% ylim([floor(min([min(subkw(:,1));min(subkw(:,2));min(subkw(:,1))])/500) ceil(max(sum(subkw, 2))/500)*500+1000])
% xlim([1 24*nptMult]);
% 
% % Plotting 24-hr flow at dist. trafo.
% figure('Renderer', 'painters', 'Position', [200 100 1100 700]);
% plot(xfmkw(:,1),'linewidth',2); hold on;
% plot(xfmkw(:,2),'linewidth',2); 
% plot(xfmkw(:,3),'linewidth',2); 
% plot(sum(xfmkw, 2),'linewidth',2);
% plot(ones(1,npts)*ratedDistTrafokw,'linewidth',2);
% hold off;
% 
% xticks(linspace(1,24*nptMult,7));
%     xticklabels({'12:00 am','04:00 am','08:00 am',...
%         '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
% grid on;
% title('Distribution transformer (634) power flow with 20% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
%     'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
% ax = gca; ax.FontSize = axisFontSize;
% str = "Sum of" + newline + "all phases";
% lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
% lgnd.FontSize = lgndFontSize; clearvars str;
% ylim([floor(min([min(xfmkw(:,1));min(subkw(:,2));min(subkw(:,1))])/500) ceil(max(sum(xfmkw, 2)))+100])
% xlim([1 24*nptMult]);
% 
% %%
% % This will contain all the losses
% % MINUS THE EV LOAD AS WELL
% % Be mindful to TAKE THE ORDER OF RUNNING
% lossesAll = [lossesAll (sum(subkw, 2) - sum(resLdkw, 2) - sum(industLdkw, 2) - sum(PevAll, 2))];
% 
% % lossesPercent = [lossesPercent sum(lossesAll(:,1))*100/sum(sum(subkw, 2))];
% 
% %% 40% EV penetration set up
% 
% close all;
% 
% extraEV = 6;
% 
% carPopKhi = (224000*0.4*0.20); % Total car sales in Pk
%                                     % times car population in Khi
%                                     % times percent penetration of EV
% workingPop = ceil(carPopKhi/100)...
%     + extraEV;                      % Car population we are working with
% 
% Vi=ones(1,numEVLds); % Initial bus reference voltage
% 
% carAtBus=readmatrix('data/evCar40pent.csv');
% 
% % Generating EV's
% [PBus, QBus] = ev_load(Vi, carAtBus);
% 
% peakTime = 12; % Kisi kaam ka nhi, hata dena isko
% carPdfMult = 1;
% 
% carPdf = normpdf(linspace(-1.1,0.5,npts),0.2,0.125) + 1*rand(1,npts); % this will multiply
%                                                                        % with the
%                                                                        % EV load to
%                                                                        % simulate gradual rise,
%                                                                        % also with some noise
% carPdf = carPdfMult*carPdf/max(carPdf);% normalizing, 
%                                            % and scaling 
%                                            % up/down with the mult
% 
% % Arrays that'll hold per-phase voltage values
% % These will be cleared for each case
% V1pu = [];
% V2pu = [];
% V3pu = [];
% 
% % power flows for substation and distribution transformer
% % These will be cleared for each case
% subkw = [];
% subkvar = [];
% xfmkw = [];
% xfmkvar = [];
% PLoss = []; 
% QLoss = [];
% PevAll = [];
% QevAll = [];
% Pev = [];
% Qev = [];
% 
% %% 40% EV penetration all plots
% 
% for i = 1:npts
%     kwLoopRes = reskw*resMult(i); kvarLoopRes = reskvar*resMult(i); % The particular power
%                                                                         % values for the i-th hour
%     kwLoopInd = industkw*industMult(i); kvarLoopInd = industkvar*industMult(i);
%     
%     if  true 
%         [Pev, Qev] = ev_load(Vi, carAtBus*carPdf(i));
%         
%         makingLoads([str_part, str_part2], [kwLoopInd(1) kwLoopRes(1:6) kwLoopInd(2:end) kwLoopRes(7:end),...
%             Pev], [kvarLoopInd(1) kvarLoopRes(1:6) kvarLoopInd(2:end) kvarLoopRes(7:end),...
%             Qev]); % make the loadss.txt file
%     else
%         makingLoads(str_part, [kwLoopInd(1) kwLoopRes(1:6) kwLoopInd(2:end) kwLoopRes(7:end)]...
%             , [kvarLoopInd(1) kvarLoopRes(1:6) kvarLoopInd(2:end) kvarLoopRes(7:end)]); % make the loadss.txt file
%     end
%     
%     PevAll = [PevAll; Pev];
%     QevAll = [QevAll; Qev];
%     
%     DailyLoadFlow; % compute the daily load flow from the files in OpenDSS
%     V1pu = [V1pu; Vpu1hr(:,1)']; % Get the per-phase voltage values for phase a
%     V2pu = [V2pu; Vpu1hr(:,2)']; % Get the per-phase voltage values for phase b
%     V3pu = [V3pu; Vpu1hr(:,3)']; % Get the per-phase voltage values for phase c
%     tempP = xlsread('MasterIEEE13_EXP_P_BYPHASE.CSV',1,'E2:J3'); % Read the powers per phase from the file
%     [tempsubkw, tempsubkvar, tempxfmkw, tempxfmkvar] = getPs(tempP); % Extract the power
%                                                                      % values for substation
%                                                                      % and distr trafo.
%     subkw = [subkw; tempsubkw]; % Store the kW power flow for substation
%     subkvar = [subkvar; tempsubkvar]; % Store the kVAR power flow for substation
%     xfmkw = [xfmkw; tempxfmkw]; % Store the kW power flow for dist. trafo.
%     xfmkvar = [xfmkvar; tempxfmkvar]; % Store the kVAR power flow for dist. trafo.
%     
%     % Putting in the phase of the load
%     Vi(1) = V1pu(i,ViSetting(1));
%     Vi(2) = V2pu(i,ViSetting(1));
%     Vi(3) = V3pu(i,ViSetting(1));
%     
%     Vi(4) = V1pu(i,ViSetting(2));
%     Vi(5) = V2pu(i,ViSetting(2));
%     Vi(6) = V3pu(i,ViSetting(2));
%     
%     Vi(7) = V1pu(i,ViSetting(3));
%     Vi(8) = V2pu(i,ViSetting(3));
%     Vi(9) = V3pu(i,ViSetting(3));
%     
%     Vi(10) = V1pu(i,ViSetting(4));
%     Vi(11) = V2pu(i,ViSetting(4));
%     Vi(12) = V3pu(i,ViSetting(4));
%     
%     Vi(13) = V1pu(i,ViSetting(5));
%     Vi(14) = V2pu(i,ViSetting(5));
%     Vi(15) = V3pu(i,ViSetting(5));
% 
%     [PLossTemp, QLossTemp] = getLinePower(readtable("MasterIEEE13_EXP_LOSSES.CSV"));
%     PLoss = [PLoss; PLossTemp];
%     QLoss = [QLoss; QLossTemp];
%     
% end
% 
% V1pu(V1pu==0) = NaN;
% V2pu(V2pu==0) = NaN;
% V3pu(V3pu==0) = NaN;
% 
% save('ev40Pent','V1pu','V2pu','V3pu',...
%     'ratedDistTrafokw','xfmkw', 'xfmkvar',...
%     'subkw','subkvar', 'ratedSubkw', 'houseLdkw','industLdkw',...
%     'houseLdkvar','industLdkvar','lossesBus','efficiencySub','efficiencyDistTrafo',...
%     'PevAll', 'QevAll');
% 
% % Load variation for all loads
% figure;
% 
% plot(sum(resLdkw, 2) + sum(industLdkw, 2), 'linewidth', plotWidth, 'Color','blue','Marker','o','MarkerFaceColor','white');  hold on
% plot(sum(resLdkw, 2) + sum(industLdkw, 2) + sum(PevAll, 2),'Color','red','Marker','o','MarkerFaceColor','white',...
%     'linewidth',plotWidth);
% grid on; xticks(linspace(1,24*nptMult,7));
% xticklabels({'12:00 am','04:00 am','08:00 am','12:00 pm','04:00 pm','08:00 pm','12:00 am'});
% title('Load variation with 40% electric vehicle penetration', 'FontSize', titleFontSize); xlabel('Time',...
%     'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
% ax = gca; ax.FontSize = axisFontSize;
% lgnd = legend(["Load with 40% EV penetration" "Normal scenario"],'Location','southoutside','orientation','horizontal');
% lgnd.FontSize = lgndFontSize;
% ylim([(floor(min(sum(resLdkw, 2) + sum(industLdkw, 2))/500)*500 - 500),...
%     (ceil(max(sum(resLdkw, 2) + sum(industLdkw, 2)+sum(PevAll, 2))/500)*500 + 500)])
% xlim([1 24*nptMult]);
% 
% % Voltages (check busNames for voltage values)
% for i =  2:length(busNames) % ignoring SOURCEBUS
%     figure('Renderer', 'painters', 'Position', [200 100 1100 700]);
%     hold on;
%     if ~isnan(V1pu(1,i))
%          plot(V1pu(:,i),'linewidth',plotWidth);
%     end
%     if ~isnan(V2pu(1,i))
%          plot(V2pu(:,i),'linewidth',plotWidth);
%     end
%     if ~isnan(V3pu(1,i))
%          plot(V3pu(:,i),'linewidth',plotWidth);
%     end
%     plot(1.05*ones(1,npts),'linewidth',plotWidth,'linestyle','--','color','black');
%     plot(0.95*ones(1,npts),'linewidth',plotWidth,'linestyle','--','color','black');
%     hold off
%     
%     xlabel('Time in hours', 'FontSize', axisFontSize);ylabel('Voltage in PU', 'FontSize', axisFontSize);
%     title(strcat("Voltage variation for bus", " ", busNames(i),", ","with 40% EV penetration (base kV_{\phi} ="," ",...
%         num2str(round(basekVs(i)/sqrt(3),3))," ","kV_{rms})"),'FontSize', titleFontSize);
%     
%     legndArr = {};
%     if ~isnan(V1pu(1,i))
%         legndArr(end+1) = {"Phase a"};
%     end
%     if ~isnan(V2pu(1,i))
%          legndArr(end+1) = {"Phase b"};
%     end
%     if ~isnan(V3pu(1,i))
%          legndArr(end+1) = {"Phase c"};
%     end
%     str1 = "1.05 p.u." + newline + "limit";
%     str2 = "0.95 p.u." + newline + "limit";
%     legndArr(end+1) = {str1};
%     legndArr(end+1) = {str2};
%     clearvars str1 str2
%     lgnd = legend(string(legndArr),'Location','southoutside','orientation','horizontal');
%     lgnd.FontSize = lgndFontSize; clear str;
%     grid on;
%     ax = gca; ax.FontSize = axisFontSize;
%     
%     xticks(linspace(1,24*nptMult,7));
%     xticklabels({'12:00 am','04:00 am','08:00 am','12:00 pm','04:00 pm','08:00 pm','12:00 am'});
%     ylim([min([V1pu(:,i)' V2pu(:,i)' V3pu(:,i)'])-0.01 1.1])
%     xlim([1 24*nptMult]);
% end
% 
% % Plotting 24-hr flow at substation
% figure('Renderer', 'painters', 'Position', [200 100 1100 700]);
% plot(subkw(:,1),'linewidth',plotWidth); hold on;
% plot(subkw(:,2),'linewidth',plotWidth); 
% plot(subkw(:,3),'linewidth',plotWidth);
% plot(sum(subkw, 2),'linewidth',plotWidth);
% plot(ones(1,npts)*ratedSubkw,'linewidth',plotWidth);
% hold off;
% 
% xticks(linspace(1,24*nptMult,7));
%     xticklabels({'12:00 am','04:00 am','08:00 am',...
%         '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
% grid on;
% title('Substation (bus 650) active power flow with 40% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
%     'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
% ax = gca; ax.FontSize = axisFontSize;
% str = "Sum of" + newline + "all phases";
% lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
% lgnd.FontSize = lgndFontSize; clearvars str;
% ylim([floor(min([min(subkw(:,1));min(subkw(:,2));min(subkw(:,1))])/500) ceil(max(sum(subkw, 2))/500)*500+1000])
% xlim([1 24*nptMult]);
% 
% % Plotting 24-hr flow at dist. trafo.
% figure('Renderer', 'painters', 'Position', [200 100 1100 700]);
% plot(xfmkw(:,1),'linewidth',plotWidth); hold on;
% plot(xfmkw(:,2),'linewidth',plotWidth); 
% plot(xfmkw(:,3),'linewidth',plotWidth); 
% plot(sum(xfmkw, 2),'linewidth',plotWidth);
% plot(ones(1,npts)*ratedDistTrafokw,'linewidth',plotWidth);
% hold off;
% 
% xticks(linspace(1,24*nptMult,7));
%     xticklabels({'12:00 am','04:00 am','08:00 am',...
%         '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
% grid on;
% title('Distribution transformer (634) power flow with 40% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
%     'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
% ax = gca; ax.FontSize = axisFontSize;
% str = "Sum of" + newline + "all phases";
% lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
% lgnd.FontSize = lgndFontSize; clearvars str;
% ylim([floor(min([min(xfmkw(:,1));min(subkw(:,2));min(subkw(:,1))])/500) ceil(max(sum(xfmkw, 2)))+100])
% xlim([1 24*nptMult]);
% 
% %%
% % This will contain all the losses
% % MINUS THE EV LOAD AS WELL
% % Be mindful to TAKE THE ORDER OF RUNNING
% lossesAll = [lossesAll (sum(subkw, 2) - sum(resLdkw, 2) - sum(industLdkw, 2) - sum(PevAll, 2))];
% 
% % lossesPercent = [lossesPercent sum(lossesAll(:,1))*100/sum(sum(subkw, 2))];
% 
% %% 90% EV penetration set up
% 
% close all;
% 
% extraEV = 439-404;
% 
% carPopKhi = (224000*0.9*0.20);  % Total car sales in Pk
%                                 % times car population in Khi
%                                 % times percent penetration of EV
% workingPop = ceil(carPopKhi/100)...
%     + extraEV;                      % Car population we are working with
% 
% Vi=ones(1,numEVLds); % Initial bus reference voltage
% 
% % Getting all EV's at the buses
% carAtBus=readmatrix('data/evCar90pent.csv');
% 
% % Generating EV's
% [PBus, QBus] = ev_load(Vi, carAtBus);
% 
% peakTime = 12; % Kisi kaam ka nhi, hata dena isko
% carPdfMult = 1;
% 
% carPdf = normpdf(linspace(-1.1,0.5,npts),0.2,0.125) + 1*rand(1,npts); % this will multiply
%                                                                        % with the
%                                                                        % EV load to
%                                                                        % simulate gradual rise,
%                                                                        % also with some noise
% carPdf = carPdfMult*carPdf/max(carPdf);% normalizing, 
%                                            % and scaling 
%                                            % up/down with the mult
% 
% % Arrays that'll hold per-phase voltage values
% % These will be cleared for each case
% V1pu = [];
% V2pu = [];
% V3pu = [];
% 
% % power flows for substation and distribution transformer
% % These will be cleared for each case
% subkw = [];
% subkvar = [];
% xfmkw = [];
% xfmkvar = [];
% PLoss = []; 
% QLoss = [];
% PevAll = [];
% QevAll = [];
% Pev = [];
% Qev = [];
% 
% %% 90% EV penetration all plots
% 
% for i = 1:npts
%     kwLoopRes = reskw*resMult(i); kvarLoopRes = reskvar*resMult(i); % The particular power
%                                                                         % values for the i-th hour
%     kwLoopInd = industkw*industMult(i); kvarLoopInd = industkvar*industMult(i);
%  
%     [Pev, Qev] = ev_load(Vi, carAtBus*carPdf(i));
% 
%     makingLoads([str_part, str_part2], [kwLoopInd(1) kwLoopRes(1:6) kwLoopInd(2:end) kwLoopRes(7:end),...
%         Pev], [kvarLoopInd(1) kvarLoopRes(1:6) kvarLoopInd(2:end) kvarLoopRes(7:end),...
%         Qev]); % make the loadss.txt file
%     
%     PevAll = [PevAll; Pev];
%     QevAll = [QevAll; Qev];
%     
%     DailyLoadFlow; % compute the daily load flow from the files in OpenDSS
%     V1pu = [V1pu; Vpu1hr(:,1)']; % Get the per-phase voltage values for phase a
%     V2pu = [V2pu; Vpu1hr(:,2)']; % Get the per-phase voltage values for phase b
%     V3pu = [V3pu; Vpu1hr(:,3)']; % Get the per-phase voltage values for phase c
%     tempP = xlsread('MasterIEEE13_EXP_P_BYPHASE.CSV',1,'E2:J3'); % Read the powers per phase from the file
%     [tempsubkw, tempsubkvar, tempxfmkw, tempxfmkvar] = getPs(tempP); % Extract the power
%                                                                      % values for substation
%                                                                      % and distr trafo.
%     subkw = [subkw; tempsubkw]; % Store the kW power flow for substation
%     subkvar = [subkvar; tempsubkvar]; % Store the kVAR power flow for substation
%     xfmkw = [xfmkw; tempxfmkw]; % Store the kW power flow for dist. trafo.
%     xfmkvar = [xfmkvar; tempxfmkvar]; % Store the kVAR power flow for dist. trafo.
%     
%     % Putting in the phase of the load
%     Vi(1) = V1pu(i,ViSetting(1));
%     Vi(2) = V2pu(i,ViSetting(1));
%     Vi(3) = V3pu(i,ViSetting(1));
%     
%     Vi(4) = V1pu(i,ViSetting(2));
%     Vi(5) = V2pu(i,ViSetting(2));
%     Vi(6) = V3pu(i,ViSetting(2));
%     
%     Vi(7) = V1pu(i,ViSetting(3));
%     Vi(8) = V2pu(i,ViSetting(3));
%     Vi(9) = V3pu(i,ViSetting(3));
%     
%     Vi(10) = V1pu(i,ViSetting(4));
%     Vi(11) = V2pu(i,ViSetting(4));
%     Vi(12) = V3pu(i,ViSetting(4));
%     
%     Vi(13) = V1pu(i,ViSetting(5));
%     Vi(14) = V2pu(i,ViSetting(5));
%     Vi(15) = V3pu(i,ViSetting(5));
% 
%     [PLossTemp, QLossTemp] = getLinePower(readtable("MasterIEEE13_EXP_LOSSES.CSV"));
%     PLoss = [PLoss; PLossTemp];
%     QLoss = [QLoss; QLossTemp];
%     
% end
% 
% V1pu(V1pu==0) = NaN;
% V2pu(V2pu==0) = NaN;
% V3pu(V3pu==0) = NaN;
% 
% save('ev90Pent','V1pu','V2pu','V3pu',...
%     'ratedDistTrafokw','xfmkw', 'xfmkvar',...
%     'subkw','subkvar', 'ratedSubkw', 'houseLdkw','industLdkw',...
%     'houseLdkvar','industLdkvar','lossesBus','efficiencySub','efficiencyDistTrafo',...
%     'PevAll', 'QevAll');
% 
% % Load variation for all loads
% figure;
% 
% plot(sum(resLdkw, 2) + sum(industLdkw, 2), 'linewidth', plotWidth, 'Color','blue','Marker','o','MarkerFaceColor','white');  hold on
% plot(sum(resLdkw, 2) + sum(industLdkw, 2) + sum(PevAll, 2),'Color','red','Marker','o','MarkerFaceColor','white',...
%     'linewidth',plotWidth);
% grid on; xticks(linspace(1,24*nptMult,7));
% xticklabels({'12:00 am','04:00 am','08:00 am','12:00 pm','04:00 pm','08:00 pm','12:00 am'});
% title('Load variation with 90% electric vehicle penetration', 'FontSize', titleFontSize); xlabel('Time',...
%     'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
% ax = gca; ax.FontSize = axisFontSize;
% lgnd = legend(["Load with 90% EV penetration" "Normal scenario" ],'Location','southoutside','orientation','horizontal');
% lgnd.FontSize = lgndFontSize;
% xlim([1 24*nptMult]);
% ylim([(floor(min(sum(resLdkw, 2) + sum(industLdkw, 2))/500)*500 - 500),...
%     (ceil(max( sum(resLdkw, 2) + sum(industLdkw, 2) + sum(PevAll, 2) )/500 )*500)])
% 
% % Voltages (check busNames for voltage values)
% for i =  2:length(busNames) % ignoring SOURCEBUS
%     figure('Renderer', 'painters', 'Position', [200 100 1100 700]);
%     hold on;
%     if ~isnan(V1pu(1,i))
%          plot(V1pu(:,i),'linewidth',plotWidth);
%     end
%     if ~isnan(V2pu(1,i))
%          plot(V2pu(:,i),'linewidth',plotWidth);
%     end
%     if ~isnan(V3pu(1,i))
%          plot(V3pu(:,i),'linewidth',plotWidth);
%     end
%     plot(1.05*ones(1,npts),'linewidth',plotWidth,'linestyle','--','color','black');
%     plot(0.95*ones(1,npts),'linewidth',plotWidth,'linestyle','--','color','black');
%     hold off
%     
%     xlabel('Time in hours', 'FontSize', axisFontSize);ylabel('Voltage in PU', 'FontSize', axisFontSize);
%     title(strcat("Voltage variation for bus", " ", busNames(i),", ","with 90% EV penetration (base kV_{\phi} ="," ",...
%         num2str(round(basekVs(i)/sqrt(3),3))," ","kV_{rms})"),'FontSize', titleFontSize);
%     
%     legndArr = {};
%     if ~isnan(V1pu(1,i))
%         legndArr(end+1) = {"Phase a"};
%     end
%     if ~isnan(V2pu(1,i))
%          legndArr(end+1) = {"Phase b"};
%     end
%     if ~isnan(V3pu(1,i))
%          legndArr(end+1) = {"Phase c"};
%     end
%     str1 = "1.05 p.u." + newline + "limit";
%     str2 = "0.95 p.u." + newline + "limit";
%     legndArr(end+1) = {str1};
%     legndArr(end+1) = {str2};
%     clearvars str1 str2
%     lgnd = legend(string(legndArr),'Location','southoutside','orientation','horizontal');
%     lgnd.FontSize = lgndFontSize; clear str;
%     grid on;
%     ax = gca; ax.FontSize = axisFontSize;
%     
%     xticks(linspace(1,24*nptMult,7));
%     xticklabels({'12:00 am','04:00 am','08:00 am','12:00 pm','04:00 pm','08:00 pm','12:00 am'});
%     ylim([min([V1pu(:,i)' V2pu(:,i)' V3pu(:,i)'])-0.01 1.1])
%     xlim([1 24*nptMult]);
% end
% 
% % Plotting 24-hr flow at substation
% figure('Renderer', 'painters', 'Position', [200 100 1100 700]);
% plot(subkw(:,1),'linewidth',plotWidth); hold on;
% plot(subkw(:,2),'linewidth',plotWidth); 
% plot(subkw(:,3),'linewidth',plotWidth);
% plot(sum(subkw, 2),'linewidth',plotWidth);
% plot(ones(1,npts)*ratedSubkw,'linewidth',plotWidth);
% hold off;
% 
% xticks(linspace(1,24*nptMult,7));
%     xticklabels({'12:00 am','04:00 am','08:00 am',...
%         '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
% grid on;
% title('Substation (bus 650) active power flow with 90% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
%     'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
% ax = gca; ax.FontSize = axisFontSize;
% str = "Sum of" + newline + "all phases";
% lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
% lgnd.FontSize = lgndFontSize; clearvars str;
% xlim([1 24*nptMult]);
% ylim([floor(min([min(subkw(:,1));min(subkw(:,2));min(subkw(:,1))])/500) ceil(max(sum(subkw, 2))/500)*500+1000])
% 
% % Plotting 24-hr flow at dist. trafo.
% figure('Renderer', 'painters', 'Position', [200 100 1100 700]);
% plot(xfmkw(:,1),'linewidth',plotWidth); hold on;
% plot(xfmkw(:,2),'linewidth',plotWidth); 
% plot(xfmkw(:,3),'linewidth',plotWidth); 
% plot(sum(xfmkw, 2),'linewidth',plotWidth);
% plot(ones(1,npts)*ratedDistTrafokw,'linewidth',plotWidth);
% hold off;
% 
% xticks(linspace(1,24*nptMult,7));
%     xticklabels({'12:00 am','04:00 am','08:00 am',...
%         '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
% grid on;
% title('Distribution transformer (634) power flow with 90% EV penetration', 'FontSize', titleFontSize); xlabel('Time',...
%     'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
% ax = gca; ax.FontSize = axisFontSize;
% str = "Sum of" + newline + "all phases";
% lgnd = legend({'Phase a','Phase b','Phase c',str,'Rated kW'},'Location','southoutside','orientation','horizontal');
% lgnd.FontSize = lgndFontSize; clearvars str;
% ylim([floor(min([min(xfmkw(:,1));min(subkw(:,2));min(subkw(:,1))])/500) ceil(max(sum(xfmkw, 2)))+100])
% xlim([1 24*nptMult]);
% 
% %%
% % This will contain all the losses
% % MINUS THE EV LOAD AS WELL
% % Be mindful to TAKE THE ORDER OF RUNNING
% lossesAll = [lossesAll (sum(subkw, 2) - sum(resLdkw, 2) - sum(industLdkw, 2) - sum(PevAll, 2))];
% 
% % lossesPercent = [lossesPercent sum(lossesAll(:,1))*100/sum(sum(subkw, 2))];
% 
% %% Losses for all scenarios
% 
% save('allLosses','lossesAll');
% 
% figure;
% plot(lossesAll,'linewidth',plotWidth);
% ylim([floor(min(min(lossesAll,[],'all'))/500) ceil(max(lossesAll,[],'all'))+30])
% xlim([1 24*nptMult]); grid on;
% xticks(linspace(1,24*nptMult,7));
% xticklabels({'12:00 am','04:00 am','08:00 am',...
%     '12:00 pm','04:00 pm','08:00 pm','12:00 am'});
% xlabel('Time'); ylabel('Varying losses (kW)');
% title('Daily losses of the system', 'FontSize', titleFontSize); xlabel('Time',...
%     'FontSize', axisFontSize); ylabel('Varying power (kW)', 'FontSize', axisFontSize);
% ax = gca; ax.FontSize = axisFontSize;
% % str = "Sum of" + newline + "all phases";
% lgnd = legend({'Normal conditions','20% EV penetration','40% EV penetration','90% EV penetration'},'Location','southoutside','orientation','horizontal');
% lgnd.FontSize = lgndFontSize; clearvars str;
% 
% %% Losses for Hammad
% 
% nptMult = 4;
% npts = 24*nptMult;
% 
% dataNormal = load('normal.mat');
% data20 = load('ev20Pent.mat');
% data40 = load('ev40Pent.mat');
% data90 = load('ev90Pent.mat');
% load('allLosses.mat');
% 
% temp1 = lossesAll(:,1)./sum(dataNormal.subkw,2)*100;
% temp2 = lossesAll(:,2)./sum(data20.subkw,2)*100;
% temp3 = lossesAll(:,3)./sum(data40.subkw,2)*100;
% temp4 = lossesAll(:,4)./sum(data90.subkw,2)*100;
% 
% maxLoss = max([temp1 temp2 temp3 temp4]);
% 
% temp1 = sum(sum(dataNormal.subkw,2))/96*24; 
% 
% plot(sum(dataNormal.subkw,2)); hold on;
% plot(sum(dataNormal.houseLdkw, 2) + sum(dataNormal.industLdkw, 2)); hold off;
% 
% sum(dataNormal.subkw,'all');
% sum(dataNormal.industLdkw ,'all') + sum(dataNormal.houseLdkw ,'all');
% 
% sum(data20.subkw,'all');
% sum(data20.industLdkw ,'all') + sum(data20.houseLdkw ,'all') + sum(data20.PevAll ,'all');
% 
% sum(data40.industLdkw ,'all') + sum(data40.houseLdkw ,'all') + sum(data40.PevAll ,'all');
% 
% sum(data90.industLdkw ,'all') + sum(data90.houseLdkw ,'all') + sum(data90.PevAll ,'all');
% %% ************************************************************************************
% 
% % Junk Code:
% % **********
% 
% % figure;
% % plot(sum(houseTemp, 2) + sum(industTemp, 2)); xticks(linspace(1,24*nptMult,7));
% %     xticklabels({'12:00 am','04:00 am','08:00 am','12:00 pm','04:00 pm','08:00 pm','12:00 am'});
% 
% % load withEVnormal;
% % area(sum(subkw, 2),'FaceColor','blue'); 
% % plot(sum(subkw, 2),'Color','blue'); %,'Marker','o','MarkerFaceColor','white');%,'-o')%,'MarkerFaceColor',[1 1 1]);
% % hold on;
% % clearvars 'V1pu' 'V2pu' 'V3pu' 'PLoss' 'QLoss' 'xfmkw' 'xfmkvar' 'subkw' 'subkvar'
% % load withoutEVnormal;
% % area(sum(subkw, 2),'FaceColor','red'); hold on;
%                                      
% % mult = [0.610897382422121	0.586227157930701	0.583430704695897...
% %     0.577909433637837	0.576945819421672	0.603618729824983	0.660501810429356...
% %     0.705171677843218	0.776959987296969	0.843342595129530	0.879595942001064...
% %     0.888725346053595	0.908307839343082	0.935718150358344	0.965077012617744...
% %     0.976324039804378	0.993509991708923	0.989814195742304	0.984771517840914...
% %     0.919952834783301	0.864804991426573	0.770438965137178	0.697671584646543...
% %     0.667986720948233];
% 
% % normal_mult = ones(1,24);
% 
% % plot(normal_mult)
% % plot(overld_mult)
% % plot(overld40_mult)
% 
% % lossesBus = ["Substation " + newline + " (650)","Distr. " + newline + " Trafo. (634)", "Line " + newline + " 650-632",...
% %     "Line " + newline + " 632-671","Line " + newline + " 671-680","Line " + newline + " 632-633","Line " + newline + " 632-645","Line " + newline + " 645-646",...
% %     "Line " + newline + " 692-675","Line " + newline + " 671-684","Line " + newline + " 684-611","Line " + newline + " 684-652","Line " + newline + " 671-692"];
% 
% % lossesBus = {"Substation (650)","Distr. Trafo. (634)", "Line 650-632",...
% %     "Line 632-671","Line 671-680","Line 632-633","Line 632-645","Line 645-646",...
% %     "Line 692-675","Line 671-684","Line 684-611","Line 684-652","Line 671-692"};
% 
% % % Plotting all active power losses
% % figure('Renderer', 'painters', 'Position', [200 100 1100 700]);
% % tempVar = [sum(PLoss,2)]';
% % bar([tempVar(1:4) tempVar(6:end-1)]) % ,'linewidth',2); hold on;
% % % xticklabels(lossesBus);
% % ax = gca; ax.FontSize = axisFontSize;
% % set(gca,'Xtick',1:11,'XtickLabel', ...
% %     [lossesBus(1:4); lossesBus(6:end-1)] , 'FontSize', axisFontSize);
% % % xlabel('Line', 'FontSize', axisFontSize);
% % ylabel('Losses (kW)', 'FontSize', axisFontSize);
% % grid on; title('Line losses for the entire network', 'FontSize', titleFontSize);
% % 
% % rotateXLabels( gca(), 30 );
% 
%%

%-----------junk code--------------
% carPopKhi = (224000*0.2*0.20);      % Total car sales in Pk
%                                     % times car population in Khi
%                                     % times percent penetration of EV
% workingPop = ceil(carPopKhi/100)...
%     + extraEV;                      % Car population we are working with
%                                     % in the network

% pf = 0.97; % power factor (assumed)

% Generating EV's
% [PBus, QBus] = ev_load(Vi, carAtBus);

% str_part2 = ["New Load.EV_634a Bus1=634.1 Phases=1 Conn=Wye  Model=1 kV=0.23 vminpu=0.3 ",...
%             'New Load.EV_634b Bus1=634.2 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
%             'New Load.EV_634c Bus1=634.3 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
%             'New Load.EV_671a Bus1=671.1 Phases=1 Conn=Wye   Model=1 kV=6.35 vminpu=0.3 ',...
%             'New Load.EV_671b Bus1=671.2 Phases=1 Conn=Wye   Model=1 kV=6.35 vminpu=0.3 ',...
%             'New Load.EV_671c Bus1=671.3 Phases=1 Conn=Wye   Model=1 kV=6.35 vminpu=0.3 ',...
%             'New Load.EV_675a Bus1=675.1 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
%             'New Load.EV_675b Bus1=675.2 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
%             'New Load.EV_675c Bus1=675.3 Phases=1 Conn=Wye   Model=1 kV=0.23 vminpu=0.3 ',...
%             'New Load.EV_692a Bus1=692.1 Phases=1 Conn=Wye   Model=1 kV=6.35 vminpu=0.3 ',...
%             'New Load.EV_692b Bus1=692.2 Phases=1 Conn=Wye   Model=1 kV=6.35 vminpu=0.3 ',...
%             'New Load.EV_692c Bus1=692.3 Phases=1 Conn=Wye   Model=1 kV=6.35 vminpu=0.3 ',...
%             'New Load.EV_680a Bus1=680.1 Phases=1 Conn=Wye   Model=1 kV=6.35 vminpu=0.3 ',...
%             'New Load.EV_680b Bus1=680.2 Phases=1 Conn=Wye   Model=1 kV=6.35 vminpu=0.3 ',...
%             'New Load.EV_680c Bus1=680.3 Phases=1 Conn=Wye   Model=1 kV=6.35 vminpu=0.3 '];


% Pavg=[11 3.7 7.2 7.4 11 11 3.7 3.7 3.7 3.7 4.6 3.2...
%     3.7 3.7 3.7 3.7 3.7 6.6 6.6 7.2 7.2 6.6 3.3 7.2 6.6...
%     6.6 7.2 7.2 3.7 11 7.2 2.8 3.7 7.2 3.7 3.7 6.6 6.6...
%     6.6 6.6 3.7 7.4 3.7 3.2 7.2 7.2 7.2 7.2 7.2 7.2 7.2...
%     3.6 3.6 7.2 3.7 3.7 22 22 22 3.3 4.6 4.6 4.6 16.5 16.5...
%     16.5 16.5 16.5 16.5 16.5 16.5 16.5 11 2.8 3.7 3.6 3.6...
%     7.2 3.6 3.6 3.6 7.2 22 3.6 3.6]; % avg. value of power in kW
%                                      % from https://www.mobilityhouse.com/int_en/
%                                      % knowledge-center/charging-time-summary

% Pev=floor(rowCopyX(PevTemp,3).*temp1./3); % The voltage dependent active power
%                                           % We are dividing by 3
%                                           % since we are loading EACH PHASE
%                                           % of an EV bus
% Qev=floor(rowCopyX(QevTemp,3).*temp2./3); % The voltage dependent reactive power

% houseLdkw = reskwDaily;
% houseLdkvar = reskvarDaily;
% industLdkw = indkwDaily;
% industLdkvar = indkvarDaily;

% efficiencySub = 0.91; % This is the efficiency from IEEE 13 bus doc
% %                                      %  We got the rated subkW by dividing
% %                                      %  total available load with
% %                                      %  efficiency
% 
% efficiencyDistTrafo = 0.984; % From the base IEEE 13 bus results

% ratedSubkw = max(sum(houseTemp, 2) + sum(industTemp, 2))/efficiencySub; % getting the rated values by dividing with efficiency
% ratedDistTrafokw = max(sum(houseTemp(:,1:3), 2))/efficiencyDistTrafo;

%----------------Junk Code---------------------

% xfmkw634(xfmkw634==ISNAN) = 0;
% xfmkvar634(xfmkvar634==ISNAN) = 0;
% 
% xfmkw675(xfmkw675==ISNAN) = 0;
% xfmkvar675(xfmkvar675==ISNAN) = 0;
% 
% xfmkw611(xfmkw611==ISNAN) = 0;
% xfmkvar611(xfmkvar611==ISNAN) = 0;
% 
% xfmkw652(xfmkw652==ISNAN) = 0;
% xfmkvar652(xfmkvar652==ISNAN) = 0;
% 
% xfmkw646(xfmkw646==ISNAN) = 0;
% xfmkvar646(xfmkvar646==ISNAN) = 0;

% plot(5500*1.1/11/sqrt(3)*ones(1,96))
% plot(5500*1.5/11/sqrt(3)*ones(1,96))
% figure; plot(subI)
% hold on
% plot(5500*1.1/11/sqrt(3)*ones(1,96))
% plot(5500*1.5/11/sqrt(3)*ones(1,96))

% houseLdkw = reskwDaily;
% houseLdkvar = reskvarDaily;
% industLdkw = indkwDaily;
% industLdkvar = indkvarDaily;

% efficiencySub = 0.91; % This is the efficiency from IEEE 13 bus doc
% %                                      %  We got the rated subkW by dividing
% %                                      %  total available load with
% %                                      %  efficiency
% 
% efficiencyDistTrafo = 0.984; % From the base IEEE 13 bus results

% ratedSubkw = max(sum(houseTemp, 2) + sum(industTemp, 2))/efficiencySub; % getting the rated values by dividing with efficiency
% ratedDistTrafokw = max(sum(houseTemp(:,1:3), 2))/efficiencyDistTrafo;

% carPdfMult = 1;
% 
% % I think the jaggedness in mult values resembles the
% % uncoordinated behavior of charging
% carpdf = normpdf(linspace(-0.8,0.5,npts),0,0.23) + 0.25*rand(1,npts); % this will multiply
%                                                                        % with the
%                                                                        % EV load to
%                                                                        % simulate gradual rise,
%                                                                        % also with some noise
% carpdf = min(carpdf)*1.25 + carPdfMult*carpdf/max(carpdf);% normalizing, 
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
% evMults = load('evMults.mat');
% carpdf = evMults.carpdf;
% bikepdf = evMults.bikepdf;
% 
% close all; plot(carpdf); hold on; % plot(circshift(carpdf, shiftA)); 
% plot(bikepdf); % plot(circshift(bikepdf, shiftB)); 
% hold off;
% xticks(linspace(1,24*nptMult,7));
%     xticklabels({'12:00 am','04:00 am','08:00 am','12:00 pm','04:00 pm','08:00 pm','12:00 am'});
%     legend("1","2");