function VpuArray1hr = getVpu1hr(csvtable)

busNames = string(csvtable.Bus);
node1 = (csvtable.Node1);
node2 = (csvtable.Node2);
% node3 = (csvtable.Node3);
V1pu1hr = (csvtable.pu1);
V2pu1hr = (csvtable.pu2);
V3pu1hr = (csvtable.pu3);

lenNode1 = length(node1);

for i = 1:lenNode1
    if node1(i) == 2
        V3pu1hr(i) = V2pu1hr(i);
        V2pu1hr(i) = V1pu1hr(i);
        V1pu1hr(i) = 0;
    else
        if node2(i) == 3
            V3pu1hr(i) = V2pu1hr(i);
            V2pu1hr(i) = 0;
        else    
            if node1(i) == 3
                V3pu1hr(i) = V1pu1hr(i);
                V2pu1hr(i) = 0;
                V1pu1hr(i) = 0;
            end
        end
    end
end

VpuArray1hr = [V1pu1hr V2pu1hr V3pu1hr];

end