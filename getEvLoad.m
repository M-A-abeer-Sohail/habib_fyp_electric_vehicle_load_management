% Total EV load
function [Pev, Qev]=getEvLoad(Vi, evAtBus, Pavg, powerFactor, busNumPh, loadMult)

nb = length(busNumPh); % Number of EV buses

nph = length(Vi); % Number of total phases that are at buses

nm = length(Pavg); % Number of EV models

p1=-0.1773; % Coefficients
p2=0.9949;
p3=0.1824;
q1=4.993;
q2=-12.91;
q3=8.917;
Vo=0.98; % nominal voltage
theta = acos(powerFactor); % phase angle

Qavg=Pavg.*tan(theta); % avg value of power in kVAR

zipTemp1 = ((p3.*(Vi/Vo).^2)+(p2.*(Vi/Vo))+p1); % here zip means the ZIP model
zipTemp2 = ((q3.*(Vi/Vo).^2)+(q2.*(Vi/Vo))+q1);

PevTemp=(evAtBus*Pavg')*loadMult; % The voltage dependent active power
QevTemp=(evAtBus*Qavg')*loadMult; % The voltage dependent reactive power
% The result of PevTemp and QevTemp, and power values with
% dimensions of nb x 1, with all power values summed at each
% bus

% A good practice is to name variables temp1, temp2 etc.,
% and then when you get their purpose, name them 
% purposefully.

% in the temp1 multiplication I need to convert both
% multiplying arrays to the same dimensions

count = 0;
Pev = [];
Qev = [];

for a = 1:nb
    perPhasePTemp = rowCopy(PevTemp(a,:), busNumPh(a)).*zipTemp1(count+1:count+busNumPh(a))'./busNumPh(a);
%     for b = busNumPh(a)
    Pev = [Pev; perPhasePTemp]; % .*(perPhasePTemp>=0)];
    perPhaseQTemp = rowCopy(QevTemp(a,:), busNumPh(a)).*zipTemp2(count+1:count+busNumPh(a))'./busNumPh(a);
%     for b = busNumPh(a)
    Qev = [Qev; perPhaseQTemp]; % .*(perPhaseQTemp>=0)];
%     end
    count = count+busNumPh(a);
end

end