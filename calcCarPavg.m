function Pavg = calcCarPavg(kWh, V, chTime)

    % Assuming charging is 80% battery capacity in kWh
    kWh = 0.9*kWh;
    Q = kWh./V;
    I = Q./chTime;
    Pavg = V.*I;

end