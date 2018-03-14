function [ trackResults ] = tracking(trackResults,settings,rawdata,loopcnt,channelNr)
%% Function is the core receiver process
% Inputs are - trackReults structure where the results are to be stored
%            - settings from where the frequency and other things are to be
%              chosen
%            - rawdata to be processed
%            - loopcnt and channelNr for which the processing is to be
%              performed
% Output - Fill the results in the structure

%% Replica code and PRN no generation for the current loop and Channel
status = trackResults(channelNr).AcqSkip(loopcnt);

if (status == 0 || status == 1 || status == 2)
ref_code = code_generator_tracking(settings,trackResults,channelNr,loopcnt);
PN_no = trackResults(channelNr).PRN;
end

%% Check Acquisition Status (0- Warm start, 2- Hot start, 1- No acquisition required  )

if (status==0)
    trackResults = Acquisition(PN_no,ref_code,trackResults,settings,rawdata,loopcnt,channelNr,status);
elseif status==2
    trackResults = Acquisition(PN_no,ref_code,trackResults,settings,rawdata,loopcnt,channelNr,status);
elseif status==1
    %% ================== Tracking Starts ========================= %%
            
% First determine the number of samples to be processed ( Corresponding to one Code or 1ms )
samples_per_code = length(ref_code); % Minimum required length for the correlation is 1 code
            
% Select the block of signal to be processed (1 ms) from recorded signal
rawSignal = rawdata((loopcnt-1)*samples_per_code+1:loopcnt*samples_per_code);

%% Carrier Wipe-off

% First select the phase and frequency to be passed to NCO for carrier
% generation
frqBin = trackResults(channelNr).carrFreq(loopcnt-1);
Phase_carr = trackResults(channelNr).phi(loopcnt-1);
% Use Carrier NCO to generate the local carrier
[carrSin,carrCos,lastPhase] = nco(frqBin,loopcnt,settings.samplingFreq,samples_per_code,Phase_carr);

% Mix the carrier with received in-coming signal (wiping-off)

            qBasebandSignal = carrSin .* rawSignal; % Quaderature
            iBasebandSignal = carrCos .* rawSignal; % In-Phase

%% Code Wipe-off and correlation computation

% Generate the replica codes (Early, Late and Prompt)

samplesPerCodeChip = round(settings.samplingFreq / settings.codeFreqBasis); % Samples in a chip
% Generate Early and Late Delays (Half Chip Shifted from the Prompt Delay)
promptDelay  = trackResults(channelNr).codephase(loopcnt-1); % Read from the previous loop
earlyDelay   = ceil(promptDelay-samplesPerCodeChip/2); % Subtract Half Chip from prompt delay
lateDelay    = ceil(promptDelay+samplesPerCodeChip/2); % Add Half Chip in prompt delay
% Generate Codes using the above Delay values
promptCode   = [ref_code(end-promptDelay+1:end) ref_code(1:end-promptDelay)];
earlyCode    = [ref_code(end-earlyDelay+1:end) ref_code(1:end-earlyDelay)];
lateCode     = [ref_code(end-lateDelay+1:end) ref_code(1:end-lateDelay)];
           
% Calculate six standard Correlation Values (Code wipe-off and accumulation)
            I_E = sum(earlyCode  .* iBasebandSignal);
            Q_E = sum(earlyCode  .* qBasebandSignal);
            I_P = sum(promptCode .* iBasebandSignal);
            Q_P = sum(promptCode .* qBasebandSignal);
            I_L = sum(lateCode   .* iBasebandSignal);
            Q_L = sum(lateCode   .* qBasebandSignal);
            
% Fill in the correlation values in results
            trackResults(channelNr).I_E(loopcnt) = I_E;
            trackResults(channelNr).Q_E(loopcnt) = Q_E;
            trackResults(channelNr).I_P(loopcnt) = I_P;
            trackResults(channelNr).Q_P(loopcnt) = Q_P;
            trackResults(channelNr).I_L(loopcnt) = I_L;
            trackResults(channelNr).Q_L(loopcnt) = Q_L;

     %% ================== Tracking Loops ========================== 
% There are two loops : Code tracking and Carrier tracking
% Carrier tracking consisits of PLL and FLL
  
%% Carrier Tracking Loop

% Phase Discriminator :

theta = atan((trackResults(channelNr).I_P(loopcnt))/(trackResults(channelNr).Q_P(loopcnt)));
trackResults(channelNr).theta(loopcnt) = theta;

% Frequency Discriminator :

% As two values of correlator outputs or phase discriminator outputs are
% required FLL loop update time is twice as that of PLL
% As FLL is used Phase values are also updated with FLL update rate
if (rem(loopcnt,2)==1) 
frqDisc = (trackResults(channelNr).theta(loopcnt)-trackResults(channelNr).theta(loopcnt-1))/(2*pi*1e-3);
%frqDisc= (atan2((trackResults(channelNr).I_P(loopcnt-1)*trackResults(channelNr).Q_P(loopcnt)-...
%trackResults(channelNr).Q_P(loopcnt-1)*trackResults(channelNr).I_P(loopcnt)),...
%(trackResults(channelNr).I_P(loopcnt-1)*trackResults(channelNr).I_P(loopcnt)+...
%trackResults(channelNr).Q_P(loopcnt-1)*trackResults(channelNr).Q_P(loopcnt))))/(2*pi);
trackResults(channelNr).freqdev(loopcnt)=frqDisc;

% Frequency Loop Filter :

trackResults(channelNr).filtfreqdev(loopcnt)=settings.k2_freq*trackResults(channelNr).filtfreqdev(loopcnt-1)+...
    settings.k1_freq*frqDisc;

% Phase Loop Filter :

trackResults(channelNr).filttheta(loopcnt)=settings.k2_phase*trackResults(channelNr).filttheta(loopcnt-2)+...
    settings.k1_phase*theta;
% Update NCO phase and frequency for next loop : 
trackResults(channelNr).phi(loopcnt) = lastPhase+trackResults(channelNr).filttheta(loopcnt);%accumulated phase for NCO
trackResults(channelNr).carrFreq(loopcnt) = frqBin+trackResults(channelNr).filtfreqdev(loopcnt);%Frequency for NCO

else % No updates at 1 ms (one code period)
trackResults(channelNr).freqdev(loopcnt)=trackResults(channelNr).freqdev(loopcnt-1);
trackResults(channelNr).filtfreqdev(loopcnt)=trackResults(channelNr).filtfreqdev(loopcnt-1);
trackResults(channelNr).phi(loopcnt) = lastPhase;
trackResults(channelNr).carrFreq(loopcnt) = frqBin;    
end

% Measure overall Doppler at every loopcnt
trackResults(channelNr).DopplerFreq(loopcnt)=trackResults(channelNr).carrFreq(loopcnt)-settings.IF;
%% Code Tracking Loop

% Calculate the Envelope values
                E = sqrt(I_E^2+Q_E^2); % Early  Envelope
                P = sqrt(I_P^2+Q_P^2); % Prompt Envelope
                L = sqrt(I_L^2+Q_L^2); % Late   Envelope
%%  E-L (Code Discriminator):
                codeError = (E - L)/(E + L);
                trackResults(channelNr).codeError(loopcnt) = codeError; 
% Check E-L and update code phase accordingly

% check if E-L curve has a zero crossing or not
if (trackResults(channelNr).codeError(loopcnt)*trackResults(channelNr).codeError(loopcnt-1)) >= 0 % No zero crossing
    
    if trackResults(channelNr).codeError(loopcnt) > 0 % Positive, shift in left
         trackResults(channelNr).codephase(loopcnt) = promptDelay-1;
    elseif trackResults(channelNr).codeError(loopcnt) < 0 % Negative, shift in right  
         trackResults(channelNr).codephase(loopcnt) = promptDelay+1;
    else % E-L=0, No shift required
         trackResults(channelNr).codephase(loopcnt) = promptDelay;
    end
    
else % zero crossing in E-L, no need to shift
   trackResults(channelNr).codephase(loopcnt) = promptDelay;
end


%% Code Phase Filter (Histogram)

if rem(loopcnt,20)==0
            hist_array=trackResults(channelNr).codephase(loopcnt-19:loopcnt); %last 20 values
            min_shift=min(hist_array); % determine max 
            max_shift=max(hist_array); % min value of code phase
            hist_length=max_shift-min_shift+1; % determine total enteries in histogram
            hist_indices = min_shift:1:max_shift; % generate indices
            hist_freq=zeros(1,hist_length); % initialize histogram
            % calculate frequency of each code phase using loops
            for i=1:20
                for j=1:hist_length
                    if hist_indices(j)==hist_array(i)
                        hist_freq(j)=hist_freq(j)+1;
                    end
                end
            end    
            [~, hist_index]=max(hist_freq); % determine maximum frequency
            trackResults(channelNr).filtcodephase(loopcnt)=hist_indices(hist_index); % update result
else
 trackResults(channelNr).filtcodephase(loopcnt)=trackResults(channelNr).filtcodephase(loopcnt-1); %pass on previous
end

            
 %% =========================== LOCK DETECTORS ===============================================

             
Lock_Detectors;             
else % any other acquisition status means satellite is not present or gone out of view
    fprintf('No Satellite Acquired for %d \n',channelNr);
    trackResults(channelNr).AcqSkip(loopcnt+1)=trackResults(channelNr).AcqSkip(loopcnt);
end

