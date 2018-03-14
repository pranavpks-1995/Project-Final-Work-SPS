%% ========== Initialize the Tracking Results ==============
%% Outputs from the correlators (In-phase):
trackResults.I_P            = zeros(1, settings.msToProcess/settings.IntegrationTime); %Prompt In-Phase
% At least 20 past values are required for SNR calculation and data demodulation 
trackResults.I_E            = zeros(1, settings.msToProcess/settings.IntegrationTime); %Early In-Phase
trackResults.I_L            = zeros(1, settings.msToProcess/settings.IntegrationTime); %Late In-Phase

% Outputs from the correlators (Quadrature-phase):
trackResults.Q_E            = zeros(1, settings.msToProcess/settings.IntegrationTime); %Prompt Quaderature
trackResults.Q_P            = zeros(1, settings.msToProcess/settings.IntegrationTime); %Early Quaderature
% At least 20 past values are required for SNR calculation and data demodulation
trackResults.Q_L            = zeros(1, settings.msToProcess/settings.IntegrationTime); %Late Quaderature

%% Carrier tracking loop parameters
trackResults.theta        = zeros(1, settings.msToProcess/settings.IntegrationTime); %Raw Phase Discriminator Output
trackResults.filttheta    = zeros(1, settings.msToProcess/settings.IntegrationTime); %Phase Loop Filter Output
trackResults.freqdev      = zeros(1, settings.msToProcess/settings.IntegrationTime); %Raw Freq Discriminator Output
trackResults.filtfreqdev  = zeros(1, settings.msToProcess/settings.IntegrationTime); %Frequency Loop Filter Output
% Carrier NCO Parameters
trackResults.phi          = zeros(1, settings.msToProcess/settings.IntegrationTime); % Accumulated Phase for NCO
trackResults.DopplerFreq  = zeros(1, settings.msToProcess/settings.IntegrationTime); % Current Doppler
trackResults.carrFreq     = ones(1, settings.msToProcess/settings.IntegrationTime)* settings.IF; % Current Carrier
%% Code Tracking Loop Parameters

% Rate at which chips of the C/A code are generated:not constant, can change due to code doppler
trackResults.codeFreq      = ones(1, settings.msToProcess/settings.IntegrationTime)* settings.codeFreqBasis;

% Tracked Code Phase (in no of samples)
trackResults.codephase     = zeros(1, settings.msToProcess/settings.IntegrationTime); 
% At least 20 past values are required for filtered code phase

% Filtered code phase(in no of samples) to be passed on to the navigation processor updated
trackResults.filtcodephase = zeros(1, settings.msToProcess/settings.IntegrationTime); % every 20 ms
trackResults.codeError     = zeros(1, settings.msToProcess/settings.IntegrationTime); %E-L sign

%% Control Lock Parameters
trackResults.LockCheck     = zeros(1, settings.msToProcess/settings.IntegrationTime); % lock count for strength
trackResults.EST_CNR       = zeros(1, settings.msToProcess/settings.IntegrationTime); % CNR Estimate (Every 20 ms)
% Hold previous 3 values for CNR
trackResults.DPrange       = zeros(1, settings.msToProcess/settings.IntegrationTime); % Every 20 ms after Lock achieved
trackResults.Data          = zeros(1, settings.msToProcess/settings.IntegrationTime);
trackResults.Data_boundary = 0;
%% Acquisition Status 

% AcqSkip 0-for warm start 1-for lock sustain 2-for lock fail (hot start required) 3-for satllite gone out of view
trackResults.AcqSkip        = zeros(1, settings.msToProcess/settings.IntegrationTime);

trackResults.AcqTh          = 0; % Decision Variable for Acquisition 
trackResults.PRN            = 0; % PRN no


%% Replicate Tracking Results Initialization for every channel
trackResults=repmat(trackResults,1,settings.numberOfChannels);

%% Cold Start, Warm Start or Hot Start
cold_start = input('Press 1 for cold start, 0 for warm start or any other for Hot Start   : ');
if cold_start==1 % Cold Start
    [trackResults, Acq_Ratio_Array] = Initial_Acquisition(trackResults,settings,rawdata);
    start_point=2;
elseif cold_start==0 % Warm Start
    start_point=1;
    for channelNr = 1:settings.numberOfChannels
        trackResults(channelNr).PRN = settings.PRN(channelNr);
    end
    
    for channelNr = 1:settings.numberOfChannels
        trackResults(channelNr).carrFreq(1:end) = settings.Init_Carr(channelNr);
    end
else % Hot Start (Update trackResults before opting for this)   
    start_point=2;
end
   
       %% ========= PRN Acquisition (warm and hot) and Tracking Loop ================== 
       
% ---------------------------------Start Processing-------------------------
    
for loopcnt=start_point:floor(settings.msToProcess/settings.IntegrationTime) % Processing for the selected time window
for channelNr = 1:settings.numberOfChannels % Processing for every channel
       loopcnt
       channelNr

   %% Computation and updation of tracking results for selected channel and loopcnt 
   trackResults = tracking(trackResults,settings,rawdata,loopcnt,channelNr);
   

end
end
%% Plot the Results when Processing is done
PlotResults