%% =========================== LOCK DETECTORS ===============================================
% The function of lock detectors is - to check lock condition,
%                                   - Update acquisition status according
%                                     to lock strength
%                                   - Demodulate data if the lock sustains
%                                   - Calculate Delta Pseudo-Range if lock
%                                     sustains
% The input of the function is the standard correlation values I_P and Q_P

%%--------------------------------- Lock Check ----------------------------   
%% SNR Criterion
if rem(loopcnt,20)==0    % lock is checked after every 20 msec
            % VSM ALGO for CNR calculation
            % Take past 20 values of I_P and Q_P
            I_P_array=trackResults(channelNr).I_P(loopcnt-19:loopcnt);
            Q_P_array=trackResults(channelNr).Q_P(loopcnt-19:loopcnt);
            % Calculate Mean and Variance of I_P^2 and Q_P^2
            Z=mean(I_P_array.^2+Q_P_array.^2);
            N=var(I_P_array.^2+Q_P_array.^2);
            % Apply VSM algo for CNR estimation
            NA=sqrt(Z^2-N); 
            est_CNR_VSM=10*log10(NA/(Z-NA)/1e-3)
            trackResults(channelNr).EST_CNR(loopcnt)=est_CNR_VSM;
%% Sign Reversal Criterion

% There may be a single data bit transition in 20 ms.
% Therefore the sign changes in the demodulated data should be at most 1 in 20 ms
% Q_P is the stream used for data in this software receiver
% therefore variable 'sign_changes_Q' should equal to 18 or more in lock
% sustain cases, but we are checking both I_P and Q_P for sign changes along 
% with  SNR and using them for determining the lock lost condition 
% and lock strength

sign_changes_I=0;
sign_changes_Q=0;
for i=1:19
sign_changes_I=sign_changes_I+sign(trackResults(channelNr).I_P(loopcnt-i)*trackResults(channelNr).I_P(loopcnt-i+1));
sign_changes_Q=sign_changes_Q+sign(trackResults(channelNr).Q_P(loopcnt-i)*trackResults(channelNr).Q_P(loopcnt-i+1));
end

%% Lock Lost        
if (((sign_changes_I<15) && (sign_changes_Q<15)) || (real(est_CNR_VSM)<30) || (imag(est_CNR_VSM)>0)) % Lock is lost

% lock may be lost if the satellite has gone out of view
% this is checked below and acq status is updated accordingly
    if (((real(est_CNR_VSM)<30) || (imag(est_CNR_VSM)>0)) && (loopcnt>60)) % check if SNR is low

        % See if SNR is decreasing gradually    
        CNR_diff1=(abs(trackResults(channelNr).EST_CNR(loopcnt))-abs(trackResults(channelNr).EST_CNR(loopcnt-20)));
        CNR_diff2=(abs(trackResults(channelNr).EST_CNR(loopcnt-20))-abs(trackResults(channelNr).EST_CNR(loopcnt-40)));
        if abs(CNR_diff1-CNR_diff2) < 2     %%% lock is lost because satellite has gone away, means no need to try again
        trackResults(channelNr).AcqSkip(loopcnt+1)=3; % Acq status is updated according to that
        trackResults(channelNr).LockCheck(loopcnt)=0;
        else % as SNR is very low - Satellite was not present,it was a false alarm  
        trackResults(channelNr).AcqSkip(loopcnt+1)=3;
        trackResults(channelNr).LockCheck(loopcnt)=0;
    end

    elseif ((((sign_changes_I<15) && (sign_changes_Q<15)) || (real(est_CNR_VSM)<30) || (imag(est_CNR_VSM)>0))...
             && loopcnt==20) % SNR is high, lock is lost for the first time due to gliches, give it a warm start
        trackResults(channelNr).AcqSkip(loopcnt+2)=0; % in the next to next loopcnt as next ms may have a transition
        trackResults(channelNr).LockCheck(loopcnt)=0;
     
    else % SNR is high, lock is lost due to gliches, give it a hot start
    trackResults(channelNr).AcqSkip(loopcnt+1)=2;
    trackResults(channelNr).LockCheck(loopcnt)=0;
    end

    trackResults(channelNr).Data(loopcnt)    = 0;
    
else % Lock is not lost, pass on previous acq status
    trackResults(channelNr).LockCheck(loopcnt)=trackResults(channelNr).LockCheck(loopcnt-1)+1;
    trackResults(channelNr).AcqSkip(loopcnt+1)=trackResults(channelNr).AcqSkip(loopcnt);
end
                    
else % wait for 20 values to accumulate, pass on previous values for lock strength and status
    trackResults(channelNr).LockCheck(loopcnt)=trackResults(channelNr).LockCheck(loopcnt-1);
    trackResults(channelNr).AcqSkip(loopcnt+1)=trackResults(channelNr).AcqSkip(loopcnt);
end

%% Delta-PseudoRange Calculation
% Delta pseudorange is the time elapsed between switching on of the
% receiver and first data bit boundary
% if the lock strength is good, check for next data transition and calculate pseudorange 
if (trackResults(channelNr).LockCheck(loopcnt)>5 && ...
sign(trackResults(channelNr).Q_P(loopcnt-1)*trackResults(channelNr).Q_P(loopcnt))==-1)

    if trackResults(channelNr).Data_boundary == 0
    trackResults(channelNr).Data_boundary = (loopcnt);
    end


    trackResults(channelNr).DPrange(loopcnt)=rem(loopcnt,20)+...
        trackResults(channelNr).filtcodephase(loopcnt)/samples_per_code;

      

    else    
    trackResults(channelNr).DPrange(loopcnt) = trackResults(channelNr).DPrange(loopcnt-1);

end

%% Data Demodulation
if (trackResults(channelNr).Data_boundary > 0 && rem((loopcnt)-trackResults(channelNr).Data_boundary,20)==0 && ...
        trackResults(channelNr).LockCheck(loopcnt)>1) % Select bit boundaries

    bit=... % fill in the value of current bit
    sign(sum(trackResults(channelNr).Q_P(loopcnt-20:loopcnt-1)));
    trackResults(channelNr).Data(loopcnt-20:loopcnt-1)=bit;

end