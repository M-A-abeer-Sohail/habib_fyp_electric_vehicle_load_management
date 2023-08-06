function [tempTrafokw, tempTrafokvar, tempLoadkw, tempLoadkvar] = getPs(tempP, LoadStart, phaseWiseVec, numTrafo)

    % We're gonna get power of all the elements for all three phases
    % Get Trafo SUB, then dist trafo 634, 675, 652, 611, and 646
    
    % The following uses some codes to 
    % arrange the load values according to phases
    % in phaseWiseVec
    % 4 refers to phase 1, delta connected load
    % 5 refers to phase 2, delta connected load
    % 6 refers to phase 3, delta connected load
    % 1 refers to phase 1, wye connected load
    % 2 refers to phase 2, wye connected load
    % 3 refers to phase 3, wye connected load
    % 0 refers to 3 phase connected load
    
    % I want to fing the row and columns where I have NaN values.
    % This gives the row and column where NaN values are present
    % The length of row and column gives the number of NaN values
    [row, col] = find(isnan(tempP));
    tempP(row, col) = 0;
    
    % Here I'm hardcoding the row from where loads start
    % but I can give these as an input to the fcn
    tempLoads = tempP(LoadStart:end,:);

    for a = 1:length(phaseWiseVec)
        if phaseWiseVec(a) == 2
            tempLoads(a, 3) = tempLoads(a, 1);
            tempLoads(a, 4) = tempLoads(a, 2);
            tempLoads(a, 1) = 0;
            tempLoads(a, 2) = 0;
        end
        if phaseWiseVec(a) == 3
            tempLoads(a, 5) = tempLoads(a, 1);
            tempLoads(a, 6) = tempLoads(a, 2);
            tempLoads(a, 1) = 0;
            tempLoads(a, 2) = 0;
        end
        if phaseWiseVec(a) == 5
            tempLoads(a, 5) = tempLoads(a, 3);
            tempLoads(a, 6) = tempLoads(a, 4);
            tempLoads(a, 3) = tempLoads(a, 1);
            tempLoads(a, 4) = tempLoads(a, 2);
            tempLoads(a, 1) = 0;
            tempLoads(a, 2) = 0;
        end
        if phaseWiseVec(a) == 6
            tempLoads(a, 5) = tempLoads(a, 1);
            tempLoads(a, 6) = tempLoads(a, 2);
            tempLoads(a, 1) = tempLoads(a, 3);
            tempLoads(a, 2) = tempLoads(a, 4);
            tempLoads(a, 3) = 0;
            tempLoads(a, 4) = 0;
        end
    end
    
    tempTrafokw = [tempP(1:numTrafo,1) tempP(1:numTrafo, 3) tempP(1:numTrafo, 5)];
    tempTrafokvar = [tempP(1:numTrafo,2) tempP(1:numTrafo, 4) tempP(1:numTrafo, 6)];
    
    tempLoadkw = [tempLoads(:,1) tempLoads(:,3) tempLoads(:,5)];
    tempLoadkvar = [tempLoads(:,2) tempLoads(:,4) tempLoads(:,6)];
    
    %------------Junk Code------------------    
%     subkw = [tempP(1,1) tempP(1,3) tempP(1,5)];
%     subkvar = [tempP(1,2) tempP(1,4) tempP(1,6)];
%     
%     xfmkw1 = [tempP(2,1) tempP(2,3) tempP(2,5)];
%     xfmkvar1 = [tempP(2,2) tempP(2,4) tempP(2,6)];
%     
%     xfmkw2 = [tempP(3,1) tempP(3,3) tempP(3,5)];
%     xfmkvar2 = [tempP(3,2) tempP(3,4) tempP(3,6)];
%     
%     xfmkw3 = [tempP(4,1) tempP(4,3) tempP(4,5)];
%     xfmkvar3 = [tempP(4,2) tempP(4,4) tempP(4,6)];
%     
%     xfmkw4 = [tempP(5,1) tempP(5,3) tempP(5,5)];
%     xfmkvar4 = [tempP(5,2) tempP(5,4) tempP(5,6)];
%     
%     xfmkw5 = [tempP(6,1) tempP(6,3) tempP(6,5)];
%     xfmkvar5 = [tempP(6,2) tempP(6,4) tempP(6,6)];
    
    %     subkw = [tempP(1,1) tempP(1,3) tempP(1,5)];
%     subkvar = [tempP(1,2) tempP(1,4) tempP(1,6)];
% 
%     xfmkw = [tempP(2,1) tempP(2,3) tempP(2,5)];
%     xfmkvar = [tempP(2,2) tempP(2,4) tempP(2,6)];
end