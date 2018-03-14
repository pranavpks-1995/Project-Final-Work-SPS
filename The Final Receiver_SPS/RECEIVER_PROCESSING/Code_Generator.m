function ref_code = Code_Generator(settings,PRN,downByFactor)
% generates reference codes according to settings, PRN and downsampling
% factor

Index=(1/(settings.samplingFreq):1/(settings.samplingFreq):settings.IntegrationTime*1e-3); % determine indices
Ind=ceil((settings.codeFreqBasis)*Index(1:end)); % Determining chip index for samples
Ind_final=rem(Ind-1,1023)+1; % Modulo 1023
    
    CAcode=generatePRN(PRN); % 1023 chips for the selected PRN
    ref_code1=CAcode(Ind_final); % sampling of chips
    ref_code=ref_code1(1:downByFactor:end); % down sampling
end