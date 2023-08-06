function Pavg = calcBikePavg(Ah, V, chTime, inkW)

    % Assuming charging is 80% battery capacity in Ah
    Ah = 0.9*Ah;
    I = Ah./chTime;
    Pavg = V.*I.*((inkW==0)*1+(inkW==1)*1/1000);

end