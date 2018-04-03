function ref_code = code_generator_tracking(settings,trackResults,channelNr,loopcnt)

%% The function gives sampled reference code considering code doppler (carrier aiding)

	Index=(1/settings.samplingFreq/2:1/settings.samplingFreq:settings.IntegrationTime*1e-3); % sampling instances  
	if loopcnt>1 % Determining the chip index for a sample
		Ind=ceil((((trackResults(channelNr).DopplerFreq(loopcnt-1))/1575.42e6+1)*settings.codeFreqBasis)*Index(1:end));%carrAid 
	
	else
		Ind=ceil((settings.codeFreqBasis)*Index(1:end));% indices of chips Initially No carrAid  
	end
	
	Ind_final=rem(Ind-1,1023)+1; % Modulo 1023 operation
	CAcode=generatePRN(trackResults(channelNr).PRN); % Code (1023 chips) generation
	ref_code=CAcode(Ind_final); % Sampling of generated Code according to generated indices

end