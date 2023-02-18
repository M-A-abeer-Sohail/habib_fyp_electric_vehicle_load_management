function [perPhaseLoad] = getDailyP_phase(allLoad, nptMult)
    % WORKS LIKE CHARM!

    % The format is as follows.
    % No of columns must be three to 
    % three phase values. Then, find
    % all phase values for each point of time.
    
    allLoadLen = size(allLoad);
    numLoad = allLoadLen(1)/(nptMult*24);
    
    perPhaseLoad = zeros(nptMult*24, 3);
    ptr = 0;
    
    for i = 1:nptMult*24
        temp = allLoad(ptr+1:ptr+numLoad,:);
        temp = sum(temp, 1);
        perPhaseLoad(i,:) = temp;
        ptr = ptr+numLoad;
    end
    
%     hello = 1;

    
    
end