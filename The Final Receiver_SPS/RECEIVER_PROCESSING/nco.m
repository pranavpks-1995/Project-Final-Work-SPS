function [ sincarr,coscarr,lastPhase ] = nco(frqBin,loopcnt,samplingFreq,samples_per_code,theta)
% Generates Sin and Cos Carriers according to supplied information 
t=1:1:samples_per_code;% Generate time vector
t = t/samplingFreq;% Indices sampling time separation
arg = frqBin*2*pi*t;% Generate arguments
arg = arg+theta;% offset the arguments by initial phase  
lastPhase=arg(end); % supply the last phase
% Generate Sin and Cos carriers
sincarr = (sin(arg));
coscarr = (cos(arg));


end

