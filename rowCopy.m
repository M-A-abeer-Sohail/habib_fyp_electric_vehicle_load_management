function copiedVec = rowCopy(toCopy, factor)
% Copies the elements of a ROW vector x amount of times

copiedVec = [];

for i = 1:factor
    copiedVec = [copiedVec; toCopy];
end

end