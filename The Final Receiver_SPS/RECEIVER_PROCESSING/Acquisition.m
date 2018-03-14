function trackResults = Acquisition(PN_no,ref_code_raw,trackResults,settings,rawdata,loopcnt,channelNr,status)
% Inputs are - Acquisition status to determine the type of acquisition
%            - Settings, loopcnt, channel no, PRN and reference code
% Updates the trackResults structure


samples_per_code = length(ref_code_raw); % No of samples taken for the processing

% The starting block of incoming data each of length equals to 1 code
signal1_raw = rawdata((loopcnt-1)*samples_per_code*settings.IntegrationTime+1:loopcnt*samples_per_code);

step_freq = 500; % The step size in frequency domain
%% Check status and start Acquisition
if status==0 % for warm start
    % find out the down sampling factor for acquisition and downsample the incoming signal  
    %downByFactor = ceil(settings.samplingFreq/settings.codeFreqBasis/12);
     downByFactor = settings.downByFactor;
    signal1 = signal1_raw(1:downByFactor:end);

    
%% code generation for reference signal    
    % Generate reference code at downsampled rate
    ref_code = ref_code_raw(1:downByFactor:end);
    samples_per_code = length(ref_code); % samples in a block after downsampling 
    
    %caCodeFreqDom = fft(fliplr(ref_code)); % FFT of replica code for Parallel Search


%% carrier generation for reference signal    
numberOfFrqBins = round(settings.acqSearchBand_warm_start * 1e3/step_freq) + 1;

%% Initialize Acquisition Results

results = zeros(numberOfFrqBins,length(ref_code));
acqRes1 = zeros(1,length(ref_code));
%--- Make the correlation for whole frequency band (for all freq. bins)
    for frqBinIndex = 1:numberOfFrqBins % Serial Frequency Search
        %frqBinIndex
        %--- Generate carrier wave frequency grid (0.5kHz step) -----------

        frqBin = trackResults(channelNr).carrFreq(loopcnt) - ...
                               (settings.acqSearchBand_warm_start/2) * 1000 + ...
                               step_freq * (frqBinIndex - 1);

        %--- Generate local sine and cosine -------------------------------
        [sinCarr,cosCarr,lastPhase] = nco(frqBin,loopcnt,settings.samplingFreq/downByFactor,length(ref_code),0);
         

        %--- "Remove carrier" from the signal -----------------------------
        I1      = sinCarr .* signal1;
        Q1      = cosCarr .* signal1;

        %% Serial Code Phase Search (comment if going for parallel search)
        
        %--- "Correlate with PN code --------------------------------------
         for i=1:length(ref_code/settings.IntegrationTime)
         acqRes1(i)=(mean([ref_code(end-i+1:end) ref_code(1:end-i)].*I1))^2+...
                (mean([ref_code(end-i+1:end) ref_code(1:end-i)].*Q1))^2;

         end
        

        %% Parallel Code Phase Search (uncomment for faster search)

       % %--- Convert the baseband signal to frequency domain --------------
        %IQfreqDom1 = fft(I1 + 1i*Q1);


      %  %--- Multiplication in the frequency domain (correlation in time
      %  %domain)
        %convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;


      %  %--- Perform inverse DFT and store correlation results ------------
        %acqRes1 = abs(ifft(convCodeIQ1)) .^ 2;

        

            results(frqBinIndex, :) = acqRes1;

        
    end
    
    
     %--- Find the correlation peak and the carrier frequency --------------
    [a, frequencyBinIndex] = max(max(results, [], 2));

    %--- Find code phase of the same correlation peak ---------------------
    [peakSize, codePhase] = max(max(results));
    
    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    samplesPerCodeChip   = round(settings.samplingFreq/downByFactor / settings.codeFreqBasis);
    excludeRangeIndex1 = codePhase - samplesPerCodeChip;
    excludeRangeIndex2 = codePhase + samplesPerCodeChip;

    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex1 < 2
        codePhaseRange = excludeRangeIndex2 : ...
                         (samples_per_code + excludeRangeIndex1);
                         
    elseif excludeRangeIndex2 >= samples_per_code
        codePhaseRange = (excludeRangeIndex2 - samples_per_code) : ...
                         excludeRangeIndex1;
    else
        codePhaseRange = [1:excludeRangeIndex1, ...
                          excludeRangeIndex2 : samples_per_code];
    end

    %--- Find the second highest correlation peak in the same freq. bin ---
    secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));
    
    if (abs(peakSize/secondPeakSize) > settings.acqThreshold) % if satellite is acquired
        %% Fine resolution frequency search ======================================        
        
        
    %--- Fill Results in structure -----  
    
        trackResults(channelNr).carrFreq(loopcnt)=trackResults(channelNr).carrFreq(loopcnt) - ...
                               (settings.acqSearchBand_warm_start/2) * 1000 + ...
                               0.5e3 * (frequencyBinIndex - 1); % Determine carr freq and update it
        % Before updating code phase, take downsampling into consideration                    
        trackResults(channelNr).codephase(loopcnt)=downByFactor*codePhase;
        trackResults(channelNr).filtcodephase(loopcnt)=downByFactor*codePhase;
        % Calculate Doppler and update
        trackResults(channelNr).DopplerFreq(loopcnt)=trackResults(channelNr).carrFreq(loopcnt)-settings.IF;
        trackResults(channelNr).AcqSkip(loopcnt+1)            = 1; % Update Acq status
        
        trackResults(channelNr).PRN            = PN_no;
        fprintf('Satellite Found %d \n ', PN_no); % Show in command
    
      %% Fine Frequency Search % uncomment if fine frequency search is needed
     reference_code = [ref_code(end-codePhase+1:end) ref_code(1:end-codePhase)];
    est_freq = trackResults(channelNr).carrFreq(loopcnt);
    fine_freq_step = 10;
    no_of_fine_frqbin = round(step_freq/fine_freq_step)+1;
    Acq_fine_result1 = zeros(1,no_of_fine_frqbin);

    for i=1:no_of_fine_frqbin
        freq = est_freq - step_freq/2 + fine_freq_step*(i-1); 
        [sinCarr,cosCarr,lastPhase] = nco(freq,loopcnt,settings.samplingFreq/downByFactor,length(ref_code),0);
        I1      = sinCarr .* signal1;
        Q1      = cosCarr .* signal1;

        Acq_fine_result1(i)=(sum(reference_code.*I1))^2+...
               (sum(reference_code.*Q1))^2;

    end
    [max1, freq_index1] = max(abs(Acq_fine_result1));

        freq_index = freq_index1;

    trackResults(channelNr).carrFreq(loopcnt)=est_freq - step_freq/2 + fine_freq_step*(freq_index-1);
    trackResults(channelNr).codephase(loopcnt)=downByFactor*codePhase;
    trackResults(channelNr).AcqSkip(loopcnt+1)            = 1;
    trackResults(channelNr).DopplerFreq(loopcnt)=trackResults(channelNr).carrFreq(loopcnt)-settings.IF; 
    
    else % satellite is not acquired
        fprintf('Satellite  Not Found %d \n ', PN_no);
        trackResults(channelNr).AcqSkip(loopcnt+1)            = 3;
    end
    trackResults(channelNr).AcqTh           = abs(peakSize/secondPeakSize);
     %% The Fine Frequency Search around estimated frequency for the ACQUISITION in case of lock fail (hot start)
elseif status==2 % Hot start 
    % Find out previous loop values
    codePhase = trackResults(channelNr).codephase(loopcnt-1);
    reference_code = [ref_code_raw(end-codePhase+1:end) ref_code_raw(1:end-codePhase)];
    est_freq = trackResults(channelNr).carrFreq(loopcnt-1);
    % choose fine frequency resolution
    fine_freq_step = 200;
    % Determine no of bins
    no_of_fine_frqbin = round(step_freq/fine_freq_step)+1;
    % Initialize hot start results
    Acq_fine_result1 = zeros(1,no_of_fine_frqbin);
   
    % start fine acquisition
    for i=1:no_of_fine_frqbin
        freq = est_freq - step_freq/2 + fine_freq_step*(i-1); 
        [sinCarr,cosCarr,lastPhase] = nco(freq,loopcnt,settings.samplingFreq,samples_per_code,0);
        I1      = sinCarr .* signal1_raw;
        Q1      = cosCarr .* signal1_raw;

        Acq_fine_result1(i)=(sum(reference_code.*I1))^2+...
               (sum(reference_code.*Q1))^2;

    end
  
    [max1, freq_index1] = max(abs(Acq_fine_result1));

        freq_index = freq_index1;

    trackResults(channelNr).carrFreq(loopcnt)=est_freq - step_freq/2 + fine_freq_step*(freq_index-1);
    trackResults(channelNr).codephase(loopcnt)=codePhase;
    trackResults(channelNr).filtcodephase(loopcnt)=codePhase;
    trackResults(channelNr).DopplerFreq(loopcnt)=trackResults(channelNr).carrFreq(loopcnt)-settings.IF;
    trackResults(channelNr).AcqSkip(loopcnt+1)            = 1;
end    
        
end

