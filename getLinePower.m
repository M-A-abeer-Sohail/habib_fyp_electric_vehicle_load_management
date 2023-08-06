function [PLoss, QLoss] = getLinePower(csvtable)

    % lineNames = string(csvtable.Element);
    Pkw = csvtable.Total_W_;
    Qkvar = csvtable.Total_var_;
    
    lenTable = 13;
    PLoss = Pkw(1:lenTable)/1000;
    QLoss = Qkvar(1:lenTable)/1000;

%           Here I took out values from POWER>CSV but they were not accurate
%           so I will use the LOSSES.CSV file 

%     PLoss = zeros(1, lenTable/2);
%     QLoss = zeros(1, lenTable/2);
%     
%     for i = 1:lenTable/2
%         PLoss(i) = abs(Pkw(2*i-1) + Pkw(2*i));
%         QLoss(i) = abs(Qkvar(2*i-1) + Qkvar(2*i));
%     end
    
end