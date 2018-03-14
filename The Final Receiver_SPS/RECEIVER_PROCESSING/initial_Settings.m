%% All the receiver parameters are listeed here.
%% Raw signal file name and other parameter ===============================
% This is the name of the data file (signal record) to be used in
% the processing 

% Select the file to be read

settings.fileName           = ...
    '/home/lenovo/Project Final Work SPS/The Final Receiver_SPS/Recorded_Signals/compactdata.bin';

%        '/home/lenovo/Project Final Work SPS/The Final Receiver_SPS/Recorded_Signals/IITJl1rxHBw_Run2.dat';
    
%        '/home/lenovo/Project Final Work SPS/The Final Receiver_SPS/Recorded_Signals/gioveAandB_short.bin';

%    '/home/lenovo/Project Final Work SPS/The Final Receiver_SPS/Recorded_Signals/compactdata.bin';

%  PRN list for warm start for the selected file

%settings.PRN=[1 3 7 11 14 19 20 22 24 28 31]; % for compactdata warm start
%settings.PRN=[22 19 14 18 32 6 11 3 28 9] % 28 9 for gioveAandB_short.bin
%start

% Initial Carrier Frequency Estimate for warm start

%settings.Init_Carr = [3566000 3559000 3566000 3564000 3563000 3560000 ...
    %3565000 3560000 3566000 3561000 3565500]; % for compactdata warm start


% Data type used to store one sample

settings.dataType           = 'ubit1'; % for compactdata.bin
%settings.dataType           = 'int8'; % for gioveAandB_short.bin
%settings.dataType           = 'int8'; % for IITJl1rxHBw_Run2.dat



% Intermediate, sampling and code frequencies

settings.IF                 = 3.563e6;      %[Hz] for compactdata.bin
%settings.IF                 = 4130400;      %[Hz] for gioveAandB_short.bin
%settings.IF                 = 6.4583475e6;      %[Hz]IITJl1rxHBw_Run2.dat


%settings.samplingFreq       = 26e6;     %[Hz]IITJl1rxHBw_Run2.dat
settings.samplingFreq       = 12e6;     %[Hz] for compactdata.bin
%settings.samplingFreq       = 16367600;     %[Hz] for gioveAandB_short.bin

settings.codeFreqBasis      = 1.023e6;      %[Hz] chip rate

% Define number of chips in a code period
settings.codeLength         = 1023;

%% Processing settings ====================================================
% Number of milliseconds to be processed

settings.msToProcess        = 500;        %[ms]

% Number of channels to be used for signal processing
settings.numberOfChannels   = 11;

% Move the starting point of processing. Can be used to start the signal
% processing at any point in the data record (e.g. for long records). fseek
% function is used to move the file read point, therefore advance is sample
% based only. 
settings.skipNumberOfBytes     = 500*26e6;
%settings.skipNumberOfBytes     =500*26e6; % for IITJl1rxHBw_Run2.dat


%% Acquisition settings ===================================================
% Band around IF to search for satellite signal. Depends on max Doppler
settings.acqSearchBand      = 20;           %[kHz] for cold start
settings.acqSearchBand_warm_start      = 10;           %[kHz] for warm start

% Threshold for the signal presence decision rule
settings.acqThreshold       = 1.9; % Considered as the threshold for the ratio of signal power to noise power

% Down sampling of the signal for acquisition
settings.downByFactor       = 1;

% The block length taken for processing
settings.IntegrationTime    = 1; % in ms, currently receiver may not work properly if it is changed

%% Tracking loops settings ================================================
% Code tracking loop parameters
% Histogram filtering is applied for code phase smoothening 

% Carrier tracking loop parameters

% FLL loop filter coefficients
settings.k1_freq=0.305;
settings.k2_freq=0.695;

% PLL loop filter coefficients
settings.k1_phase=1;
settings.k2_phase=0;



