%% After Processing is over, plot the results for different satellites

No_of_channels = settings.numberOfChannels;

%% Acquisition Results Plot

if cold_start==1
Acqed_Ratio_Array = zeros(1,length(Acq_Ratio_Array));
for i=1:length(Acq_Ratio_Array)
    if Acq_Ratio_Array(i) > settings.acqThreshold
     Acqed_Ratio_Array(i) = Acq_Ratio_Array(i);
    end
end
    figure(500), bar(Acq_Ratio_Array,'r'), hold on, bar(Acqed_Ratio_Array,'g')
end


for channelNr=1:No_of_channels
    %% create figure for all the channels
    figure(channelNr);
    clf(channelNr);
    set(channelNr, 'Name', ['Channel ', num2str(channelNr), ...
                                 ' (PRN ', ...
                                 num2str(trackResults(channelNr).PRN), ...
                                 ') Acquisition Ratio',num2str(trackResults(channelNr).AcqTh)]);
                             
   %% Arrange subplots
        % Row 1
        handles(1, 1) = subplot(3, 3, 1);
        handles(1, 2) = subplot(3, 3, 2);
        handles(1, 3) = subplot(3, 3, 3);
        % Row 2
        handles(2, 1) = subplot(3, 3, 4);
        handles(2, 2) = subplot(3, 3, 5);
        handles(2, 3) = subplot(3, 3, 6);
        % Row 3
        handles(3, 1) = subplot(3, 3, 7);
        handles(3, 2) = subplot(3, 3, 8);
        handles(3, 3) = subplot(3, 3, 9);
    %% Plot the Tracking Results
    % plot time(in seconds) on x-axis
      timevector = (1:settings.msToProcess)/1000;
    %% Plot code phase values
     plot  (handles(1, 1), timevector, ...
                              (trackResults(channelNr).filtcodephase)/(settings.samplingFreq*1e-3));

        grid  (handles(1, 1));
        title (handles(1, 1), 'Filtered Code Phase');
        xlabel(handles(1, 1), 'Time (s)');
        axis  (handles(1, 1), 'tight');
     %% Plot Doppler Frequency Values
     plot  (handles(1, 2), timevector, ...
                              (trackResults(channelNr).DopplerFreq));

        grid  (handles(1, 2));
        title (handles(1, 2), 'Doppler Frequency');
        xlabel(handles(1, 2), 'Time (s)');
        axis  (handles(1, 2), 'tight');
        %% Plot filtered frequency deviation
        plot  (handles(1, 3), timevector, ...
                              (trackResults(channelNr).filtfreqdev));

        grid  (handles(1, 3));
        title (handles(1, 3), 'FLL Filter Output');
        xlabel(handles(1, 3), 'Time (s)');
        axis  (handles(1, 3), 'tight');
        %% Plot filtered phase deviation
        plot  (handles(2, 1), timevector, ...
                              (trackResults(channelNr).filttheta));

        grid  (handles(2, 1));
        title (handles(2, 1), 'PLL Filter Output');
        xlabel(handles(2, 1), 'Time (s)');
        axis  (handles(2, 1), 'tight');
        %% Plot I_P and Q_P values
        plot  (handles(2, 2), timevector, ...
                              (trackResults(channelNr).I_P),'r');
                          
        hold on
        
        plot  (handles(2, 2), timevector, ...
                              (trackResults(channelNr).Q_P));
                          
        grid  (handles(2, 2));
        title (handles(2, 2), 'Prompt Corr O/P');
        xlabel(handles(2, 2), 'Time (s)');
        axis  (handles(2, 2), 'tight');
        %% Plot Data values
        plot  (handles(2, 3), timevector, ...
                              (trackResults(channelNr).Data));

        grid  (handles(2, 3));
        title (handles(2, 3), 'Demodulated Data');
        xlabel(handles(2, 3), 'Time (s)');
        axis  (handles(2, 3), 'tight');
        %% Lock Fail Check
        scatter  (handles(3, 1), timevector, ...
                              (trackResults(channelNr).AcqSkip(1:end-1)));

        grid  (handles(3, 1));
        title (handles(3, 1), 'Lock Fails (at 2)');
        xlabel(handles(3, 1), 'Time (s)');
        axis  (handles(3, 1), 'tight');
        %% Lock Fail Check
        plot  (handles(3, 2), timevector, ...
                              (trackResults(channelNr).DPrange));
        
        grid  (handles(3, 2));
        title (handles(3, 2), 'Delta Pseudo-range');
        xlabel(handles(3, 2), 'Time (s)');
        axis  (handles(3, 2), [timevector(1) timevector(end) trackResults(channelNr).DPrange(end)-0.001 trackResults(channelNr).DPrange(end)+0.001]);
        
        %% Lock Fail Check
        plot  (handles(3, 3), timevector(20:20:end), ...
                              (trackResults(channelNr).EST_CNR(20:20:end)));

        grid  (handles(3, 3));
        title (handles(3, 3), 'Estimated CNR');
        xlabel(handles(3, 3), 'Time (s)');
        axis  (handles(3, 3), 'tight');
end  