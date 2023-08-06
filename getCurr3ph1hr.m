function Curr3ph1hr = getCurr3ph1hr(tempTable, startRow, endRow)
    % tempTable should be of 6 columns
    
    Curr3ph1hr = [tempTable(startRow:endRow,1) tempTable(startRow:endRow,3) tempTable(startRow:endRow,5)];
end