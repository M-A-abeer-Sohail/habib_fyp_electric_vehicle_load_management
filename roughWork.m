nb = 3;

nph = 5;

nm = 4;

busNumPh = [3,1,1];
Vi = ones(1,nph);

evAtBus = ones(nb, nm);

Pavg = ones(1,nm);

pf = 0.9;

getEvLoad(Vi, evAtBus, Pavg, pf, busNumPh)
%------------------------------------------------------------
% temp = getDailyP_phase(industLdkw, nptMult);
% plot(sum(temp, 2));

% temp = ones(96*8,3);
% temp2 = getDailyP_phase(temp, nptMult);
% reshape(.',1,[]) % Just a syntax I guess
% reshape(ans,[8,3])

% temp = xlsread('MasterIEEE13_EXP_P_BYPHASE.CSV',1,'E2:J34');
% 
% % This gives the row and column where NaN values are present
% % The length of row and column gives the number of NaN values
% [row, col] = find(isnan(temp));
% temp(row, col) = 0;
% 
% % Here I'm hardcding the row from where loads start
% % but I can give these as an input to the fcn
% tempLoads = temp(16:end,:);
% 
% % The following uses some codes to 
% % arrange the load values according to phases
% % 4 refers to phase 1, delta connected load
% % 5 refers to phase 2, delta connected load
% % 6 refers to phase 3, delta connected load
% % 1 refers to phase 1, wye connected load
% % 2 refers to phase 2, wye connected load
% % 3 refers to phase 3, wye connected load
% % 0 refers to 3 phase connected load
% phaseWise = [4 5 6 0 0 1 2 3 1 2 3 1 2 1 3 1 2 3];
% 
% for a = 1:length(phaseWise)
%     if phaseWise(a) == 2
%         tempLoads(a, 3) = tempLoads(a, 1);
%         tempLoads(a, 4) = tempLoads(a, 2);
%         tempLoads(a, 1) = 0;
%         tempLoads(a, 2) = 0;
%     end
%     if phaseWise(a) == 3
%         tempLoads(a, 5) = tempLoads(a, 1);
%         tempLoads(a, 6) = tempLoads(a, 2);
%         tempLoads(a, 1) = 0;
%         tempLoads(a, 2) = 0;
%     end
%     if phaseWise(a) == 5
%         tempLoads(a, 5) = tempLoads(a, 3);
%         tempLoads(a, 6) = tempLoads(a, 4);
%         tempLoads(a, 3) = tempLoads(a, 1);
%         tempLoads(a, 4) = tempLoads(a, 2);
%         tempLoads(a, 1) = 0;
%         tempLoads(a, 2) = 0;
%     end
%     if phaseWise(a) == 6
%         tempLoads(a, 5) = tempLoads(a, 1);
%         tempLoads(a, 6) = tempLoads(a, 2);
%         tempLoads(a, 1) = tempLoads(a, 3);
%         tempLoads(a, 2) = tempLoads(a, 4);
%         tempLoads(a, 3) = 0;
%         tempLoads(a, 4) = 0;
%     end
% end

