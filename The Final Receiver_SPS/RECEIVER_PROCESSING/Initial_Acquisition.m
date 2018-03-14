function [trackResults, Acq_Ratio_Array] = Initial_Acquisition(trackResults,settings,rawdata)
% Function is used for a cold start case
% Inputs are rawdata and settings
% Takes rawdata, perform acquisition and updates trackResults 
% find out the no of satellites visible and acquisition parameters for them


    No_of_visible_sats=0; % initialize the variables
    
    samples_per_code = (settings.samplingFreq/settings.codeFreqBasis)*1.023e3;% No of samples taken for the processing
    
    % The starting blocks of incoming data each of length equals to 1 code
    signal1_raw = rawdata(1:samples_per_code);

    %% Downsampling of the incoming signal for acquisition
    %downByFactor = ceil(settings.samplingFreq/settings.codeFreqBasis/12);
    downByFactor = settings.downByFactor;
    samples_per_code = (settings.samplingFreq/downByFactor/settings.codeFreqBasis)*1.023e3;%new no of samples per block
    
    step_freq = 500;% The step size in frequency domain
    % Downsample 
    signal1 = signal1_raw(1:downByFactor:end);

   
    
    PRN=0; % Initially
    PRN_Max = 32; % Max possible PRN for GPS change it for IRNSS
    Acq_Ratio_Array = zeros(1,PRN_Max); %For storing acq ratio of all the satellites 
    % Perform Acquisition for the no of channels available
    
    for channelNr=1:settings.numberOfChannels
        
    AcqTh = 0; % Initailization
    
    % Search Every PRN in increasing order and go to next channel whenever
    % one is acquired
    while (AcqTh < settings.acqThreshold && PRN < PRN_Max)
        
    PRN = PRN+1; % Increment in PRNs
    
    % Ref code for selected PRN
    ref_code = Code_Generator(settings,PRN,downByFactor);
    
    %caCodeFreqDom = fft(fliplr(ref_code));% FFT for parallel search
    numberOfFrqBins = round(settings.acqSearchBand * 1e3/step_freq) + 1;% No of bins in frequency
    
    %Initialize acq results
    results = zeros(numberOfFrqBins,length(ref_code));
    acqRes1 = zeros(1,length(ref_code));
%--- Make the correlation for whole frequency band (for all freq. bins)
    for frqBinIndex = 1:numberOfFrqBins % Serial Frequency Search
        %frqBinIndex
        %--- Generate carrier wave frequency grid (0.5kHz step) -----------
        frqBin = settings.IF - ...
                               (settings.acqSearchBand/2) * 1000 + ...
                               step_freq * (frqBinIndex - 1);

        %--- Generate local sine and cosine -------------------------------
 [sinCarr,cosCarr,~] = nco(frqBin,1,settings.samplingFreq/downByFactor,length(ref_code),0);
         

        %--- "Remove carrier" from the signal -----------------------------
        I1      = sinCarr .* signal1;
        Q1      = cosCarr .* signal1;

        
        %% Serial Code Phase Search (comment if choosing parallel code phase search)
       % %--- "Correlate with PN code"--------------------------------------
        for i=1:length(ref_code/settings.IntegrationTime)
        acqRes1(i)=(mean([ref_code(end-i+1:end) ref_code(1:end-i)].*I1))^2+...
               (mean([ref_code(end-i+1:end) ref_code(1:end-i)].*Q1))^2;

        end
        

        
       %% Parallel Code Phase Search (Uncomment for faster processing)
       
       % %--- Convert the baseband signal to frequency domain --------------
        %IQfreqDom1 = fft(I1 + 1i*Q1);


       % %--- Multiplication in the frequency domain (correlation in time
       % %domain)
        %convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;


       % %--- Perform inverse DFT and store correlation results ------------
        %acqRes1 = abs(ifft(convCodeIQ1)) .^ 2;

        

      %  % Fill the results for a frequency bin in Grid
            results(frqBinIndex, :) = acqRes1;

        
    end
    
    
     %--- Find the correlation peak and the carrier frequency --------------
    [~, frequencyBinIndex] = max(max(results, [], 2));

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
    
    AcqTh = (abs(peakSize/secondPeakSize));
    
    Acq_Ratio_Array(PRN) = AcqTh;
    
    if AcqTh > settings.acqThreshold
        
        
    %--- Fill Results in structure -----    
        trackResults(channelNr).carrFreq(1)=settings.IF - ...
                               (settings.acqSearchBand/2) * 1000 + ...
                               0.5e3 * (frequencyBinIndex - 1);
        trackResults(channelNr).codephase(1)=downByFactor*codePhase;
        trackResults(channelNr).filtcodephase(1)=downByFactor*codePhase;
        trackResults(channelNr).promptDelay(1)=downByFactor*codePhase;
        trackResults(channelNr).DopplerFreq(1)=trackResults(channelNr).carrFreq(1)-settings.IF;
        trackResults(channelNr).AcqSkip(2)            = 1;
        
        trackResults(channelNr).Status            = 'T';
        trackResults(channelNr).PRN            = PRN;
         trackResults(channelNr).AcqTh            = AcqTh;
        fprintf('Satellite Found %d \n ', PRN);
    
      %% Fine Frequency Search 
     reference_code = [ref_code(end-codePhase+1:end) ref_code(1:end-codePhase)];
    est_freq = trackResults(channelNr).carrFreq(1);
    fine_freq_step = 10;
    no_of_fine_frqbin = round(step_freq/fine_freq_step)+1;
    Acq_fine_result1 = zeros(1,no_of_fine_frqbin);

    for i=1:no_of_fine_frqbin
        freq = est_freq - step_freq/2 + fine_freq_step*(i-1);
[sinCarr,cosCarr,lastPhase] = nco(freq,1,settings.samplingFreq/downByFactor,length(ref_code),0);
        I1      = sinCarr .* signal1;
        Q1      = cosCarr .* signal1;

        Acq_fine_result1(i)=(sum(reference_code.*I1))^2+...
               (sum(reference_code.*Q1))^2;

    end
    [max1, freq_index1] = max(abs(Acq_fine_result1));

        freq_index = freq_index1;

    trackResults(channelNr).carrFreq(1)=est_freq - step_freq/2 + fine_freq_step*(freq_index-1);
    
    trackResults(channelNr).AcqSkip(2)            = 1;
    trackResults(channelNr).AcqTh            = AcqTh;
    trackResults(channelNr).DopplerFreq(1)=trackResults(channelNr).carrFreq(1)-settings.IF;
    No_of_visible_sats = No_of_visible_sats+1;
    end
    end
    end
    % If sufficient no of satellites not visible, disable the remaining channels
    
    if No_of_visible_sats < settings.numberOfChannels
        for channelNr=No_of_visible_sats+1:settings.numberOfChannels
            trackResults(channelNr).AcqSkip(2)            = 3;
        end
    end
    
end