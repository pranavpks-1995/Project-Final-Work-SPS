%% Clean up the environment first =========================================
close all; clc; 

format ('compact');
format ('long', 'g');
%% Startup Printing ======================================================= 
fprintf(['\n',...
    'Welcome to:  THE Software GNSS Receiver\n\n', ...
    'A GNSS software project by:\n\n', ...
    '              Wireless Communication Laboratory\n', ...
    '             Department of Electrical Engineering\n',...
    '            Indian Institute of Technology Jodhpur\n']);
fprintf('                -------------------------------\n\n');
%% Initialize Settings =========================================
%----Read from initial settings file and store in workspace   
initial_Settings

%% Recorded Data File Initialization =========================================================
%----Read the recorded data file-----
disp ('Starting processing...');

[fid, message] = fopen(settings.fileName, 'rb','b');

%If success, then process the data
if (fid > 0)
        
    % Set the starting point for the processing in the recorded data file 
    fseek(fid, settings.skipNumberOfBytes, 'bof');
    samplesPerCode = round(settings.samplingFreq / ...
                           (settings.codeFreqBasis / settings.codeLength));
    % Read data for processing 
    rawdata = fread(fid, settings.msToProcess * samplesPerCode , settings.dataType)';
        %-----Modified to read new data file
        rawdata = rawdata * 2 - 1; %To convert values to +/- 1 for compactdata.bin file only
else
    % Error while opening the data file.
    error('Unable to read file %s: %s.', settings.fileName, message);
end 
% Close the data file
  fclose(fid);
  %% Generate plot of raw data and ask if ready to start processing =========
  probe_Data = input('Press 1 if you want to probe the data \n 0 for directly going to processing ');
  
  if probe_Data==1
      
fprintf('Probing data (%s)...\n', settings.fileName)
probeData(settings);

disp('  Raw IF data plotted ')
disp('  (run setSettings or change settings in "initial_Settings.m" to reconfigure)')
disp(' ');

  end
  
Start = input('Enter "1" to initiate the processing or "0" to exit : ');
if (Start == 1)
    disp(' ');
    %start things rolling...
    Processing
end
