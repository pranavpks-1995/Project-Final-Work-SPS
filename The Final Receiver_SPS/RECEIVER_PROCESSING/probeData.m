function probeData(settings)
% Function plots the PSD of rawdata. 
% It is used for checking if the data record is fine or problematic.
%% Generate plot of raw data ==============================================
[fid, message] = fopen(settings.fileName, 'rb');


if (fid > 0)
    % Move the starting point of processing. Can be used to start the
    % signal processing at any point in the data record (e.g. for long
    % records).
    fseek(fid, settings.skipNumberOfBytes, 'bof');    
    
    % Find number of samples per spreading code
    samplesPerCode = round(settings.samplingFreq / ...
                           (settings.codeFreqBasis / settings.codeLength));
                      
    % Read 10ms of signal
    [data, count] = fread(fid, [1, 10*samplesPerCode], settings.dataType);
    
    
    fclose(fid);
    
    if (count < 10*samplesPerCode)
        % The file is too short
        error('Could not read enough data from the data file.');
    end
    
 %% After reading the file plot PSD
    figure(51)
    
    pwelch(data-mean(data), 16384, 1024, 2048, settings.samplingFreq/1e6)
else
    %=== Error while opening the data file ================================
    error('Unable to read file %s: %s.', fileNameStr, message);
end % if (fid > 0)
