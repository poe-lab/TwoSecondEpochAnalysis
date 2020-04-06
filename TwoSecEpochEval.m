Fs= 1024; % This is a constant for this program.
    
% EEG bandwidths:
D_lo = 0.4; % Specify low end of Delta band
D_hi = 4.9; % Specify high end of Delta band
T_lo = 5; % Specify low end of Theta band
T_hi = 9;   % Specify high end of Theta band
S_lo = 10; % Specify low end of Sigma band
S_hi = 14;  % Specify high end of Sigma band

% EEG filter set:
EEG_Fc = 30; % The EEG low-pass cut-off frequency. Default is 30 Hz.
EEG_highpass_enable = 0; % Set to '0' to turn off high-pass filter, '1' to turn on.
EEG_HP_Fc = 1; % The EEG high-pass cut-off frequency. Default is 1 Hz.
EEG_Notch_enable = 0; % Set to '0' to turn off notch filter, '1' to turn on.

% EMG filter set:
EMG_Fc = 30; % The EMG high-pass cut-off frequency. Default is 30 Hz.
EMG_Notch_enable=0; % Set to '0' to turn off notch filter, '1' to turn on.
EMG_lowpass_enable=0;  % Set to '0' to turn off low-pass filter, '1' to turn on.

for i = 1:numOfEpochs
    st_pt=EPOCH_StartPoint;
    end_pt=st_pt+EPOCHSIZE-1; 
    Vj=double(EMG_SAMPLES(st_pt:end_pt));
    absVj=abs(Vj).^2;                       % absolute square
    P_emg(INDEX)=sum(absVj)/length(absVj);  % sum of all squared Vj's
    
    % This is calculating power in frequency domain
    fft_in=double(EEG_SAMPLES(st_pt:end_pt));
    windowsize =length(fft_in);
    if df < windowsize
        df = windowsize;
    end
    
%     h = spectrum.welch('Hann', ones(windowsize,1), 0);  % Form: h = spectrum.welch('Hann',window,100*noverlap/window);
%     hpsd = psd(h, double(fft_in), 'NFFT', df, 'Fs', Fs);    %Form: hpsd = psd(h,x,'NFFT',nfft,'Fs',Fs);
%     Pxx2 = hpsd.Data;
%     F2 = hpsd.Frequencies;
    % SPECTRUM is obsolete. Replaced with above 4 lines that does the equivalent.
    [Pxx2,F2]=spectrum(double(fft_in),df,0,ones(windowsize,1),Fs);
    % ******  [P,F] = SPECTRUM(X,NFFT,NOVERLAP,WINDOW,Fs)
    
    %For the EEG signal
    index_delta=[];index_theta=[];index_sigma=[]; index_beta=[];
    index_delta=find(F2(1)+D_lo< F2 & F2 < F2(1)+D_hi);      % Default delta band 0.4 -4 Hz
    index_theta=find(F2(1)+T_lo< F2 & F2 < F2(1)+T_hi);    % Default theta band 5-9 Hz
    index_sigma=find(F2(1)+S_lo< F2 & F2 < F2(1)+S_hi);     % Default sigma band 10-14 Hz
    index_beta =find(F2(1)+B_lo< F2 & F2 < F2(1)+B_hi);     % Default Beta band 15-20 Hz
    
    P_delta(INDEX)=sum(Pxx2(index_delta))/df *2;  
    P_theta(INDEX)=sum(Pxx2(index_theta))/df *2;
    P_sigma(INDEX)=sum(Pxx2(index_sigma))/df *2;
    P_beta(INDEX)=sum(Pxx2(index_beta))/df *2;
    
    st_power(INDEX)=abs(P_sigma(INDEX)*P_theta(INDEX));   % Used to indicate waking
    dt_ratio(INDEX)=abs(P_delta(INDEX)/P_theta(INDEX));   
    
    Statetime(INDEX)=EEG_TIMESTAMPS(EPOCH_StartPoint);
    EPOCH_StartPoint=st_pt+EPOCHSIZE;        % Keep this at the end of this function
    INDEX=INDEX+1;
    warning('off','MATLAB:divideByZero');
    
    
end

% New Parameter on next line:
averageSTpower = mean(st_power);
StdDevforAverageSTpower = std(st_power);

