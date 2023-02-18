% clear all;

clc;
DSSObj = actxserver('OpenDSSEngine.DSS');
if ~DSSObj.Start(0),
disp('Unable to start the OpenDSS Engine');
return
end

% delete('MasterIEEE13_EXP_METERS.CSV')
delete('MasterIEEE13_EXP_P_BYPHASE.CSV');
delete('MasterIEEE13_EXP_VOLTAGES.CSV');
delete('MasterIEEE13_EXP_POWERS.CSV');
delete('MasterIEEE13_EXP_LOSSES.CSV');

DSSText = DSSObj.Text; % Used for all text interfacing from matlab to opendss
DSSCircuit = DSSObj.ActiveCircuit; % active circuit

tempTextCommand=['Compile (',fileparts(which('DailyLoadFlow.m')),'\MasterIEEE13.dss)']; % Get's the 
                                                                                        % path to the 
                                                                                        % file, wherever 
                                                                                        % it is.
DSSText.Command = tempTextCommand;% Path where Master and its associated files are stored.
DSSText.Command='batchedit load..* Vmin=0.3'; % Set Vmin so that load model property will remain same
DSSTransformers=DSSCircuit.Transformers;
nt=24*1;
TimeArray=1:nt;

DSSText.Command='set mode=snapshot ! daily stepsize=1h number=1';
DSSText.Command='set hour=0'; % Start at second 0 of hour
DSSText.Command='set BaseFrequency=50'; % Start at second 0 of hour

    DSSText.Command = 'Solve';
    DSSText.Command = 'Export p_byphase';
    DSSText.Command = 'Export voltages';
    DSSText.Command = 'Export powers';
    DSSText.Command = 'Export losses';
    DSSText.Command = 'Export currents';
    
    temp = readtable('MasterIEEE13_EXP_VOLTAGES.csv');
    
    Vpu1hr = getVpu1hr(temp);
%% Plotting voltages

% busNames = 1:12;
% busNames = string(busNames);
% 
% sizeVol = size(V1pu);
% for i = 1:sizeVol(2)
%     figure('Renderer', 'painters', 'Position', [5 5 1400 800]);
% plot(V1pu(:,i)); hold on;
% plot(V2pu(:,i));
% plot(V3pu(:,i)); hold off
% xlabel('Time in hrs');ylabel('Voltage in PU');
% title({'Phase voltages at bus', busNames(i)});
% xticks(1:24);
%     xticklabels({'01','02','03','04','05','06','07','08','09'...
%         '10','11','12','13','14','15','16','17','18',...
%         '19','20','21','22','23','00'});
%     legend('Phase a','Phase b','Phase c');
% % saveas(gcf,'Varying Active Loads for bus 5, 6 and 8.png');
% end

%% Messy Version
% %%**************
% 
% clear all;
% clc;
% DSSObj = actxserver('OpenDSSEngine.DSS');
% if ~DSSObj.Start(0),
% disp('Unable to start the OpenDSS Engine');
% return
% end
% 
% delete('MasterIEEE13_EXP_METERS.CSV')
% delete('MasterIEEE13_EXP_P_BYPHASE.CSV');
% delete('MasterIEEE13_EXP_VOLTAGES.CSV');
% 
% mult = [0.610897382422121	0.586227157930701	0.583430704695897...
%     0.577909433637837	0.576945819421672	0.603618729824983	0.660501810429356...
%     0.705171677843218	0.776959987296969	0.843342595129530	0.879595942001064...
%     0.888725346053595	0.908307839343082	0.935718150358344	0.965077012617744...
%     0.976324039804378	0.993509991708923	0.989814195742304	0.984771517840914...
%     0.919952834783301	0.864804991426573	0.770438965137178	0.697671584646543	0.667986720948233];
% 
% DSSText = DSSObj.Text; % Used for all text interfacing from matlab to opendss
% DSSCircuit = DSSObj.ActiveCircuit; % active circuit
% DSSText.Command='Compile (C:\Users\Abeer\OneDrive - Habib University\V2G Potential FYP\FYP codes & data\SatsangiMatlab\MasterIEEE13.dss)';% Path where Master and its associated files are stored.
% DSSText.Command='batchedit load..* Vmin=0.8'; % Set Vmin so that load model property will remain same
% DSSTransformers=DSSCircuit.Transformers;
% %DSSText.Command='Batchedit regcontrol..* Enabled=no'; % uncomment for tap change as per user's choice
% % DSSText.Command='batchedit load..* daily=1'; % Loadshape
% % DSSText.Command='New EnergyMeter.Main Line.650632 1';% Energy meter
% nt=24*1;
% TimeArray=1:nt;
% 
% %% Uncomment following for Tap chang as per user choice (manually)
% % Xtap=[15	7	6	6	6	7	8	9	11	12	13	14	14	14	14	14	14	14	14	14	14	14	13	12
% % 10	5	4	4	4	4	5	6	7	8	8	9	9	9	9	9	9	9	9	9	9	9	9	9
% % 15	6	5	5	5	6	7	9	10	12	13	13	14	14	14	14	14	14	14	14	14	14	13	12 ];
% % Reg1Tap=Xtap(1,:)-5;
% % Reg2Tap=Xtap(2,:)-5;
% % Reg3Tap=Xtap(3,:)-5;
% % Vreg1=1+0.00625*Reg1Tap;
% % Vreg2=1+0.00625*Reg2Tap;
% % Vreg3=1+0.00625*Reg3Tap;
% DSSText.Command='set mode=snapshot ! daily stepsize=1h number=1';
% % DSSText.Command='set mode= daily stepsize=1h number=1';
% DSSText.Command='set hour=0'; % Start at second 0 of hour
% % for i=1:nt
% % DSSText.Command='get hour';
% %    hour=DSSText.Result;
% 
%     %% Uncomment following for change in Tap Positions of Regulator as per user choice
%     % DSSText.command = ['Transformer.Reg1.Tap=',num2str(Vreg1(i))];
%     % DSSText.command = ['Transformer.Reg2.Tap=',num2str(Vreg2(i))];
%     % DSSText.command = ['Transformer.Reg3.Tap=',num2str(Vreg3(i))];
%     DSSText.Command='Solve';
% %     SystemLosses(i,:)=(DSSCircuit.Losses)/1000; % Will Give you Distribution System Losses in kWs and kVArs
% %     %% Line Losses
% %     LineLosses(i,:)=DSSCircuit.Linelosses;
% %     %% Transformer Losses
% % 
% %     TranLosses(i,:)=SystemLosses(i,:)-LineLosses(i,:);
% % 
% %     %% Voltage Magnitude in p.u. for 24-hours can be obtained in this way
% %     % V1pu(i,:)=DSSCircuit.AllNodeVmagPUByPhase(1);
% %     % V2pu(i,:)=DSSCircuit.AllNodeVmagPUByPhase(2);
% %     % V3pu(i,:)=DSSCircuit.AllNodeVmagPUByPhase(3);
% %     DSSText.Command = '? Transformer.Reg1.Taps';
% %     Reg1=str2num(DSSText.Result); 
% %     Vreg1S(i,:)=Reg1(2);% Secondary winding voltage of Reg1 in 24-hr
% %     DSSText.Command = '? Transformer.Reg2.Taps';
% %     Reg2=str2num(DSSText.Result);
% %     Vreg2S(i,:)=Reg2(2);% Secondary winding voltage of Reg2 in 24-hr
% %     DSSText.Command = '? Transformer.Reg3.Taps';
% %     Reg3=str2num(DSSText.Result); 
% %     Vreg3S(i,:)=Reg3(2);% Secondary winding voltage of Reg3 in 24-hr
% %    DSSText.Command='Export Meter'; % A MasterIEEE13_EXP_METERS.CSV file will be saved in same path
%     DSSText.Command='Export p_byphase';
%     DSSText.Command = 'Export voltages';
%     
%     temp = readtable('MasterIEEE13_EXP_VOLTAGES.csv');
%     
%     Vpu1hr = getVpu1hr(temp);
% %     V1pu(i,:) = Vpu1hr(:,1)';
% %     V2pu(i,:) = Vpu1hr(:,2)';
% %     V3pu(i,:) = Vpu1hr(:,3)';
% % end
% 
% % DSSText.Command = 'Export voltages';
% % EM=csvread('MasterIEEE13_EXP_METERS.CSV',1,4);
% % SubkWh=EM(:,1);
% % SubkVArh=EM(:,2);
% % SubkW24=[SubkWh(TimeArray(1)); SubkWh(TimeArray(2):TimeArray(end))-SubkWh(TimeArray(1):TimeArray(nt-1))];
% % SubkVAr24=[SubkVArh(TimeArray(1)); SubkVArh(TimeArray(2):TimeArray(end))-SubkVArh(TimeArray(1):TimeArray(nt-1))];
% % SubkVA24=abs(SubkW24+sqrt(-1)*SubkVAr24);
% % SubstationkWkVArandkVA24=[SubkW24 SubkVAr24 SubkVA24];
% 
% % Here SubkW24 SubkVAr24 and SubkVA24 are substation kWs, kVArs, and kVAs for 24 hours period
% % clearvars -except temp SystemLosses LineLosses SubkW24 SubkVAr24 SubkVA24 V1pu V2pu V3pu Vreg1S Vreg2S Vreg3S
% 
% %% Plotting voltages
% 
% % busNames = 1:12;
% % busNames = string(busNames);
% % 
% % sizeVol = size(V1pu);
% % for i = 1:sizeVol(2)
% %     figure('Renderer', 'painters', 'Position', [5 5 1400 800]);
% % plot(V1pu(:,i)); hold on;
% % plot(V2pu(:,i));
% % plot(V3pu(:,i)); hold off
% % xlabel('Time in hrs');ylabel('Voltage in PU');
% % title({'Phase voltages at bus', busNames(i)});
% % xticks(1:24);
% %     xticklabels({'01','02','03','04','05','06','07','08','09'...
% %         '10','11','12','13','14','15','16','17','18',...
% %         '19','20','21','22','23','00'});
% %     legend('Phase a','Phase b','Phase c');
% % % saveas(gcf,'Varying Active Loads for bus 5, 6 and 8.png');
% % end
% 
% %%
% % clc;
% % temp = readtable('MasterIEEE13_EXP_VOLTAGES.CSV');
% % busNames = string(temp.Bus);
% % node1 = (temp.Node1);
% % node2 = (temp.Node2);
% % node3 = (temp.Node3);
% % V1pu1hr = (temp.pu1);
% % V2pu1hr = (temp.pu2);
% % V3pu1hr = (temp.pu3);
% % 
% % lenNode1 = length(node1);
% % 
% % for i = 1:lenNode1
% %     if node1(i) == 2
% %         V3pu1hr(i) = V2pu1hr(i);
% %         V2pu1hr(i) = V1pu1hr(i);
% %         V1pu1hr(i) = 0;
% %     else
% %         if node2(i) == 3
% %             V3pu1hr(i) = V2pu1hr(i);
% %             V2pu1hr(i) = 0;
% %         else    
% %             if node1(i) == 3
% %                 V3pu1hr(i) = V1pu1hr(i);
% %                 V2pu1hr(i) = 0;
% %                 V1pu1hr(i) = 0;
% %             end
% %         end
% %     end
% % end
% 
% % temp2 = [V1pu1hr V2pu1hr V3pu1hr];

