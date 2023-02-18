% str1 = "New Load.671 Bus1=671.1.2.3  Phases=3 Conn=Delta Model=1 kV=4.16   ";
% str2 = 1155; str3 = 660;
% 
% str_dss = strcat(str1,"kw = ",num2str(str2), "   kvar = ",num2str(str3));

% str_part = ["New Load.671 Bus1=671.1.2.3  Phases=3 Conn=Delta Model=1 kV=4.16   ",...
%     "New Load.634a Bus1=634.1     Phases=1 Conn=Wye   Model=1 kV=0.277"  ];
% kw = [1155 160];
% kvar = [660 110];

function makingLoads(str_part, kw, kvar)

    fid = fopen('Loadss.txt','w');
    
    str_dss = "! LOAD DEFINITIONS FOR IEEE-13 BUS RADIAL DISTRIBUTION SYSTEM";
    fprintf(fid,'%s \n',str_dss);

    maxCount = length(kw);
    
%     str_dss = str_part(1);
%         fprintf(fid,'%s \n',str_dss);
    
    for a = 1:maxCount
%         if a == 8 || a == 9 || a == 10 || a == 1
%             str_dss = str_part(a);
%         else
            str_dss = strcat(str_part(a),"kw = ",num2str(floor(kw(a))), "   kvar = ",num2str(floor(kvar(a))));
%         end
        fprintf(fid,'%s \n',str_dss);
    end

    fclose(fid);
end