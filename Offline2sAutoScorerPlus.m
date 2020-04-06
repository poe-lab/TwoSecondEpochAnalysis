function Offline2sAutoScorerPlus
%--------------------------------------------------------------------------
% DEFINE CONSTANTS

Fs= 1000; % This is a constant for this program.
    
% EEG bandwidths:
D_lo = 0.4; % Specify low end of Delta band
D_hi = 4; % Specify high end of Delta band
T_lo = 4.1; % Specify low end of Theta band
T_hi = 9;   % Specify high end of Theta band
S_lo = 9.1; % Specify low end of Sigma band
S_hi = 15;  % Specify high end of Sigma band
B_lo = 15; 
B_hi = 20;

% EEG filter set:
EEG_Fc = 30; % The EEG low-pass cut-off frequency. Default is 30 Hz.
EEG_highpass_enable = 0; % Set to '0' to turn off high-pass filter, '1' to turn on.
EEG_HP_Fc = 1; % The EEG high-pass cut-off frequency. Default is 1 Hz.
EEG_Notch_enable = 0; % Set to '0' to turn off notch filter, '1' to turn on.

% EMG filter set:
EMG_Fc = 30; % The EMG high-pass cut-off frequency. Default is 30 Hz.
EMG_Notch_enable=0; % Set to '0' to turn off notch filter, '1' to turn on.
EMG_lowpass_enable=0;  % Set to '0' to turn off low-pass filter, '1' to turn on.

% -------------------------------------------------------------------------
% SELECT DATA FILES:

% Select EMG file
working_dir=pwd;
current_dir='C:\SleepData\DataFiles';
cd(current_dir);
[fileName, pathName] = uigetfile({'*.ncs','Neuralynx CSC File (*.ncs)'},...
    'Select the EMG data file');
emgFile= fullfile(pathName, fileName);
cd(working_dir);
clear fileName pathName

% Select EEG file
working_dir=pwd;
current_dir='C:\SleepData\DataFiles';
cd(current_dir);
[fileName, pathName] = uigetfile({'*.ncs','Neuralynx CSC File (*.ncs)'},...
    'Select the EEG data file');
eegFile= fullfile(pathName, fileName);
cd(working_dir);
clear fileName pathName

% Select the time stamp file 
working_dir=pwd;
current_dir='C:\SleepData\Timestampfiles';
cd(current_dir);
[fileName, pathName] = uigetfile('*.xls', 'Select the time stamp file');
timeStampFile = fullfile(pathName, fileName);
cd(working_dir);
clear fileName pathName

% Select the manually scored file that is in 10s epochs and convert to 2s epochs 
[manualScored2sEpochs, scoredTS, fileLabel] = User10to2sScoringConverter;
numOf10sEpochs = size(scoredTS,2);
numOf2sEpochs = size(manualScored2sEpochs,2);


%Write header for PSD file before entering loop to calculate PSd values for
%all 2s epochs
c = clock;
dt = datestr(c,'mmddyy-HHMM');
resultExcelFileName = ['PsdAndAS2sEpochsBasedOn' fileLabel '_' dt '.xls'];
warning off MATLAB:xlswrite:AddSheet
sheetName = 'Files used';
xlswrite(resultExcelFileName,cellstr('EMG'), sheetName, 'A1');
xlswrite(resultExcelFileName,cellstr('EEG'), sheetName, 'A2');
xlswrite(resultExcelFileName,cellstr('Time Stamps'), sheetName, 'A3');
xlswrite(resultExcelFileName,cellstr('User Scored'), sheetName, 'A4');
xlswrite(resultExcelFileName,cellstr(emgFile), sheetName, 'B1');
xlswrite(resultExcelFileName,cellstr(eegFile), sheetName, 'B2');
xlswrite(resultExcelFileName,cellstr(timeStampFile), sheetName, 'B3');
xlswrite(resultExcelFileName,cellstr(fileLabel), sheetName, 'B4');

        

% Read in the timestampsfilebutton file.
[tbounds] = xlsread(timeStampFile);
lbound = tbounds(1:end,1);  ubound = tbounds(1:end,2);
exactLowIndx = tbounds(1:end,3); exactHiIndx = tbounds(1:end,4);

Sum_P_sigma=0; Length_P_sigma=0; Squaresum_P_sigma=0;
Mean_sigma=[]; Std_dev_sigma=[];
scored2sEpochTs = [];
INDEX = [];
warning('off', 'signal:spectrum:obsoleteFunction');
for boundIndex=1:length(ubound)
    lowertimestamp=lbound(boundIndex);
    uppertimestamp=ubound(boundIndex);
    exactLow = exactLowIndx(boundIndex);
    exactHi = exactHiIndx(boundIndex);
    % Thresholds is made '[]' because its only initialized when we
    % autoscore it.. Its a structure having all threshold values
    Thresholds=[];
    nsamp = []; % Number of samples per bin of time

    %Import EMG Data:
    [Timestamps,SF,Samples] = Nlx2MatCSC(emgFile,[1 0 1 0 1],0,4,[lowertimestamp uppertimestamp]);
    samples = double(Samples(:)');
    clear Samples
    SampFreq1 = SF(1);
    if boundIndex == 1
        DS = (1:1:10);
        DSampSF = SampFreq1./DS;
        indSampfactor = find(DSampSF >= 1000);
        Fs = round(DSampSF(indSampfactor(end)));
        sampfactor = DS(indSampfactor(end));
         msgbox({['Orginal Sampling Rate:  ' num2str(SampFreq1) 'Hz'];...
            ['Down-Sampled Sampling Rate:  ' num2str(Fs) 'Hz']; ['Sampling Factor:  ' num2str(sampfactor) '']});
    end
    % Precise time stamps should be calculated here:
    [EMG_TIMESTAMPS,EMG_SAMPLES] = generate_timestamps_from_Ncsfiles(Timestamps, samples, exactLow, exactHi, nsamp);
    %Set up EMG filters:
    %  High pass filter for EMG signals
    [Bhigh,Ahigh]=ellip(7,1,60, EMG_Fc/(SampFreq1/2),'high');  % Default setting implements high pass filter with 30hz cutoff
    physInput = 1;  %Needed to select proper error box in HeaderADBit.
    ADBit2uV = HeaderADBit(emgFile, physInput);    %Calls a function to extract the AD Bit Value.
    EMG_SAMPLES = EMG_SAMPLES * ADBit2uV;   %Convert EMG amplitude of signal from AD Bits to microvolts.
    filtered_samples = filter(Bhigh,Ahigh,EMG_SAMPLES);
    % OPTIONAL lowpass filter for EMG signals
    if EMG_lowpass_enable>0
        [EMG_Blow,EMG_Alow] = ellip(7,1,60, EMG_LP_Fc/(SampFreq1/2));   % Default is OFF
        filtered_samples = filter(EMG_Blow,EMG_Alow, filtered_samples);
    end
    % Optional 60Hz Notch Filter
    if EMG_Notch_enable > 0
        woB = 60/(SampFreq1/2);
        [B_EMG_Notch,A_EMG_Notch] =  iirnotch(woB, woB/35);   % Default is OFF
        filtered_samples = filter(B_EMG_Notch,A_EMG_Notch, filtered_samples);
    end
    EMG_TIMESTAMPS = EMG_TIMESTAMPS(1:sampfactor:end);
    EMG_SAMPLES = filtered_samples(1:sampfactor:end);
    clear physInput ADBit2uV

    clear filtered_samples samples

    %  ******    EEG FILE extraction   *********
    [Timestamps,SF,Samples]=Nlx2MatCSC(eegFile,[1 0 1 0 1],0,4,[lowertimestamp uppertimestamp]);
    samples=double(Samples(:)');
    clear Samples
    SampFreq2=SF(1);
    Fs=round(SampFreq2/sampfactor);
    [EEG_TIMESTAMPS,EEG_SAMPLES] = generate_timestamps_from_Ncsfiles(Timestamps,samples, exactLow, exactHi, nsamp);
    physInput = 2;  %Needed to select proper error box in HeaderADBit.
    ADBit2uV = HeaderADBit(eegFile, physInput);    %Calls a function to extract the AD Bit Value.
    EEG_SAMPLES = EEG_SAMPLES * ADBit2uV;   %Convert EEG amplitude of signal from AD Bits to microvolts.
    %  Low pass filter for EEG signals
    [Blow,Alow]=ellip(7,1,60, EEG_Fc/(SampFreq2/2));           % Default setting implements low pass filter with 30hz cutoff
    filtered_samples=filter(Blow,Alow,EEG_SAMPLES);
    %  OPTIONAL highpass filter for EEG signals
    if EEG_highpass_enable>0
        [EEG_Bhi,EEG_Ahi] = ellip(7,1,60, EEG_HP_Fc/(SampFreq2/2),'high');   % Default is OFF
        filtered_samples = filter(EEG_Bhi,EEG_Ahi, filtered_samples);
    end
    %  OPTIONAL 60Hz Notch filter for EEG signals
    if EEG_Notch_enable > 0
        wo = 60/(SampFreq2/2);
        [B_EEG_Notch,A_EEG_Notch] =  iirnotch(wo, wo/35);   % Default is OFF
        filtered_samples = filter(B_EEG_Notch,A_EEG_Notch, filtered_samples);
    end
    EEG_TIMESTAMPS = EEG_TIMESTAMPS(1:sampfactor:end);
    EEG_SAMPLES = filtered_samples(1:sampfactor:end);
    clear physInput ADBit2uV

    clear filtered_samples samples SF
    %*******************************************************************
    scored10sIndexFor2HrBin = find(scoredTS >= EMG_TIMESTAMPS(1)-1 & scoredTS <= EMG_TIMESTAMPS(end)-5);
    num10sEpochsFor2HrBin = length(scored10sIndexFor2HrBin);
    
    for i = 1:num10sEpochsFor2HrBin-1
        if isequal(i,720)
            bob=1;
        end
        startTime = scoredTS(scored10sIndexFor2HrBin(i));
        index = find(EMG_TIMESTAMPS > startTime - 0.01 & EMG_TIMESTAMPS < startTime + 0.01);
        diff= EMG_TIMESTAMPS(index(1):index(end)) - startTime;
        [~,ind]=min(abs(diff));
        indexStartTimeInEmgTS = index(ind);
        for j = 1:5
            startPoint = indexStartTimeInEmgTS + (j-1)*2*Fs;
            stopPoint = startPoint + 2*Fs - 1;
            ts2sEpoch = EMG_TIMESTAMPS(startPoint);
            scored2sEpochTs = [scored2sEpochTs; ts2sEpoch];
            INDEX = size(scored2sEpochTs,1);
            emg2sData = EMG_SAMPLES(startPoint:stopPoint);
            eeg2sData = EEG_SAMPLES(startPoint:stopPoint);
            
            % EMG power calculation for 2s epoch
            Vj=double(emg2sData);
            absVj=abs(Vj).^2;                       % absolute square
            powerEMG(INDEX,1) =sum(absVj)/length(absVj);  % sum of all squared Vj's
            
            % This is calculating EEG power in frequency domain
            fft_in=double(eeg2sData);
            windowsize =length(fft_in);
            df = 2*Fs;
            if df < windowsize
                df = windowsize;
            end

            % SPECTRUM is obsolete. Replaced with above 4 lines that does the equivalent.
            [Pxx2,F2]=spectrum(fft_in,df,0,ones(windowsize,1),Fs);
            % ******  [P,F] = SPECTRUM(X,NFFT,NOVERLAP,WINDOW,Fs)

            %For the EEG signal
            index_delta=[];index_theta=[];index_sigma=[]; index_beta=[];
            index_delta=find(F2(1)+D_lo< F2 & F2 < F2(1)+D_hi);      % Default delta band 0.4 -4 Hz
            index_theta=find(F2(1)+T_lo< F2 & F2 < F2(1)+T_hi);    % Default theta band 5-9 Hz
            index_sigma=find(F2(1)+S_lo< F2 & F2 < F2(1)+S_hi);     % Default sigma band 10-14 Hz
            index_beta =find(F2(1)+B_lo< F2 & F2 < F2(1)+B_hi);     % Default Beta band 15-20 Hz

            deltaPower(INDEX,1)=sum(Pxx2(index_delta))/df *2;  
            thetaPower(INDEX,1)=sum(Pxx2(index_theta))/df *2;
            sigmaPower(INDEX,1)=sum(Pxx2(index_sigma))/df *2;
            betaPower(INDEX,1)=sum(Pxx2(index_beta))/df *2;
            clear index_delta index_theta index_sigma index_beta
            st_power(INDEX,1)=abs(sigmaPower(INDEX,1).*thetaPower(INDEX,1));   % Used to indicate waking
            dt_ratio(INDEX,1)=abs(deltaPower(INDEX)./thetaPower(INDEX));   
            warning('off','MATLAB:divideByZero');
        end      
    end
    clear EMG_SAMPLES EMG_TIMESTAMPS EEG_SAMPLES EEG_TIMESTAMPS  
end
numOf2sEpochs = INDEX;
%Write PSD information to a new sheet in the Excel file.
sheetName = 'PSDs-2sEpochs';
xlswrite(resultExcelFileName,cellstr('Timestamp'), sheetName, 'A1');
xlswrite(resultExcelFileName,cellstr('MS Stage'), sheetName, 'B1');
xlswrite(resultExcelFileName,cellstr('EMG Power'), sheetName, 'C1');
xlswrite(resultExcelFileName,cellstr('Delta Power'), sheetName, 'D1');
xlswrite(resultExcelFileName,cellstr('Theta Power'), sheetName, 'E1');
xlswrite(resultExcelFileName,cellstr('Sigma Power'), sheetName, 'F1');
xlswrite(resultExcelFileName,cellstr('Beta Power'), sheetName, 'G1');

xlswrite(resultExcelFileName,[scored2sEpochTs manualScored2sEpochs powerEMG deltaPower thetaPower sigmaPower betaPower],sheetName, 'A2');


% Calculates the Mean and Std Deviation for sigma:
meanSigmaPower = mean(sigmaPower);
stdDevSigmaPower = std(sigmaPower);

thresholdStatus = 0;
def = {'','','','','','',''};
while thresholdStatus < 7 %#ok<NODEF>
    thresholdStatus = 0;
    % Dialog box to enter the power thresholds for auto-scoring:
    prompt = {'Enter EMG Threshold:','Enter SxT Threshold:','Enter D/T Threshold:',...
        'Enter Delta Threshold:','Enter Theta Threshold:',...
        'Enter Sigma Threshold:','Enter Sigma Standard Deviation:'};
    dlg_title = 'Auto-scoring Thresholds';
    lineNo = 1;

    answer = inputdlg(prompt,dlg_title,lineNo,def);
    def = answer';
    for i = 1:7 %Get all of the entered threshold values and make sure they are numbers.
        if isequal(answer{i,1}, '') %Check to see if empty threshold box
        elseif isnan(str2double(answer{i,1}))   %Check to see if not a number
        else
            Thresholds(i) = str2double(answer{i,1});    %#ok<AGROW> %It is a number, so add to threshold vector
            thresholdStatus = thresholdStatus + 1;
        end
    end
    clear answer prompt dlgTitle lineNo
end

unhookedThreshold = 0.001;
% choiceSaveThresholds = questdlg('Save threshold values?', '',...
%          'Yes','No','Yes');
% 
% switch choiceSaveThresholds
%     case 'Yes'
%         %Request user input to name time stamp file:
%         prompt = {'File name for threshold values:'};
%         def = {'SubjectNumber'};
%         dlgTitle = 'Save Threshold Settings';
%         lineNo = 1;
%         answer = inputdlg(prompt,dlgTitle,lineNo,def);
%         filename = char(answer(1,:));
%         dateString = date;
%         thresholdFilename = strcat('C:\SleepData\', filename, dateString, '.xls');
% 
%         thresholdArray = {'Threshold Settings for 2 second epochs', dateString;...
%             'EMG:', Thresholds(1); 'Sigma*Theta:', Thresholds(2); 'Delta/Theta:', Thresholds(3);...
%             'Delta:',Thresholds(4); 'Theta:', Thresholds(5); 'Sigma:', Thresholds(6);...
%             'Sigma Std Dev:', Thresholds(7); 'Unhooked:', unhookedThreshold};
% 
%         xlswrite(thresholdFilename,thresholdArray)
%         clear prompt def dlgTitle lineNo answer filename thresholdFilename thresholdArray
%     case 'No'
%         % Continue to experiment phase.
% end

autoScoredEpochs = [];
pastWake = [1 1 0];
%--------------------------------------------------------------------------
% FULL AUTOSCORER LOGIC 
for i= 1:numOf2sEpochs
    
    tempState = 2;  % Set the default state as QS.

    % Score epochs either as REM or QS depending on their DT ratio and EMG
    if dt_ratio(i) < Thresholds(3) && powerEMG(i) < Thresholds(1)
        tempState = 3;  % Change the state to REM.
    end

    % WAKING decision point
    if powerEMG(i) > Thresholds(1) && st_power(i) < Thresholds(2)
    % Now from all those epochs which can be termed as QW/ AW based on
    % STthresh
        tempState = 1;  % Change the state to AW.
    end

    % Check if the threshold is given & then look for the UNhooked epochs  
    % The logic in this section is not performed if no threshold is entered.
    if isempty(unhookedThreshold)==0  
        if powerEMG(i) < 1 && st_power(i) < unhookedThreshold
            tempState = 5;
        end
    end

    % For absolute detection of AW states from states with very high EMG as
    % well as high Sigma * Theta value which are scored as QS before this
    if isequal(tempState,2)
        if powerEMG(i) > Thresholds(1) && st_power(i) > Thresholds(2)
            tempState = 1;  % Change to AW
        else
            % Check if they have low D/T power and should be REM
            if dt_ratio(i) < Thresholds(3) && powerEMG(i) < (1.1*Thresholds(1));
                tempState = 3;  % Change to REM
            end
        end
    end

    % Check if there are any TR state within QS states
    if isequal(tempState,2)
%                 if st_power > (5*averageSTpower)    % + 10*StdDevforAverageSTpower));
%                     tempState = 6;  % Change to TR
%                 end
        %Original published logic:
        if sigmaPower(i) > (meanSigmaPower + 2*stdDevSigmaPower)  % This line is in the analysis of 1 s epochs.
            % The multiplier is 3SD in 10s epochs.
            tempState = 6;  % Change to TR
        end

    end             

% To see if there is any QW within the AW or QS states by checking
% their sigma, delta and theta levels. They should be below the
% threshold set by the user
    if isequal(tempState,1) || (isequal(tempState,2) && deltaPower(i) < Thresholds(4))
        if thetaPower(i) < Thresholds(5) && sigmaPower(i) < Thresholds(6)
            if powerEMG(i) < 2.5*Thresholds(1)
                tempState = 4;  % Change to QW.  This logic may need to be reviewed.
            end
        end
    end

    % For absolute REM detection with the REM states
    if isequal(tempState,3)
        %Change from REM to Qiuet Wake if the last three epochs were AW (=1) or
        %QW (=4)
        if isequal(pastWake, [1 1 1])
            tempState = 4;
        end
    end

% This will not be used unless we want to add new logic to make IW 
% detection a determinant for another state.
    if isequal(tempState,4)
        if thetaPower(i) > Thresholds(5)
            % Change to Intermediate Wake (IW) for correction.
            tempState = 8;
        end
    end

    %Save the auto-scored state
    autoScoredEpochs = [autoScoredEpochs; tempState];
    if isequal(tempState, 1) || isequal(tempState, 4) || isequal(tempState, 8)
        isWakeState = 1;
    else
        isWakeState = 0;
    end
    pastWake = [pastWake(2:3) isWakeState];
end     
%Write Auto-Scored epoch information to a new sheet in the Excel file.
sheetName = 'AS-2sEpochs';
xlswrite(resultExcelFileName,cellstr('Timestamp'), sheetName, 'A1');
xlswrite(resultExcelFileName,cellstr('AS Stage'), sheetName, 'B1');

xlswrite(resultExcelFileName,[scored2sEpochTs autoScoredEpochs],sheetName, 'A2');


%--------------------------------------------------------------------------
%                           COMPARISON MATRIX

agree = 0;

% Find total percent agreement
for k = 1:numOf2sEpochs
    if isequal(manualScored2sEpochs(k), autoScoredEpochs(k))
        agree = agree + 1;
    end
end
percentAgree = agree/numOf2sEpochs;

t=1;
numberMismatch = zeros(6,7);

for i = 1:8
    if i == 7
    else
        agree = 0;
        index = find(manualScored2sEpochs(:) == i);
        p = length(index);
        numP(i) = p;
        manualState = manualScored2sEpochs(index);
        if isempty(manualState) == 0
            autoState = autoScoredEpochs(index);

            for k = 1:p
                if isequal(manualState(k), autoState(k))
                    agree = agree + 1;
                else
                     switch autoState(k)
                        case 1
                            numberMismatch(t,1) =  numberMismatch(t,1) + 1;
                        case 2
                            numberMismatch(t,2) =  numberMismatch(t,2) + 1;
                        case 3
                            numberMismatch(t,3) =  numberMismatch(t,3) + 1;
                        case 4
                            numberMismatch(t,4) =  numberMismatch(t,4) + 1;
                        case 5
                            numberMismatch(t,5) =  numberMismatch(t,5) + 1;
                        case 6
                            numberMismatch(t,6) =  numberMismatch(t,6) + 1;
                        case 8
                            numberMismatch(t,7) =  numberMismatch(t,7) + 1;
                     end
                end
            end
        numberMismatch(t,t) = agree;
        end
        t = t + 1;
    end
end

for i = 1:6
    manualTotals(i) = sum(numberMismatch(i,:)); % User
    %autoTotals(i) = sum(numberMismatch(:,i)); % Auto-Scorer
    manualVsAuto(i) = numberMismatch(i,i)/manualTotals(i);
    %autoVsManual(i) = numberMismatch(i,i)/autoTotals(i);
end
%autoTotals(7) = sum(numberMismatch(:,7));  %Intermediate Waking total scored by Auto-Scorer.

%Write performance information to a new sheet in the Excel file.
sheetName = 'AS Performance';

xlswrite(resultExcelFileName, cellstr('Total % Agreement'), sheetName, 'A1');
xlswrite(resultExcelFileName, percentAgree, sheetName, 'B1');

xlswrite(resultExcelFileName,cellstr('AW'), sheetName, 'B3');
xlswrite(resultExcelFileName,cellstr('QS'), sheetName, 'C3');
xlswrite(resultExcelFileName,cellstr('RE'), sheetName, 'D3');
xlswrite(resultExcelFileName,cellstr('QW'), sheetName, 'E3');
xlswrite(resultExcelFileName,cellstr('UH'), sheetName, 'F3');
xlswrite(resultExcelFileName,cellstr('TR'), sheetName, 'G3');
xlswrite(resultExcelFileName,cellstr('IW'), sheetName, 'H3');
xlswrite(resultExcelFileName,cellstr('Manual Total'), sheetName, 'I3');
xlswrite(resultExcelFileName,cellstr('AS % Agree'), sheetName, 'J3');

xlswrite(resultExcelFileName,cellstr('AW'), sheetName, 'A4');
xlswrite(resultExcelFileName,cellstr('QS'), sheetName, 'A5');
xlswrite(resultExcelFileName,cellstr('RE'), sheetName, 'A6');
xlswrite(resultExcelFileName,cellstr('QW'), sheetName, 'A7');
xlswrite(resultExcelFileName,cellstr('UH'), sheetName, 'A8');
xlswrite(resultExcelFileName,cellstr('TR'), sheetName, 'A9');
xlswrite(resultExcelFileName,cellstr('IW'), sheetName, 'A10');

xlswrite(resultExcelFileName,[numberMismatch(1:6,:) manualTotals' manualVsAuto'],sheetName, 'B4');

%Calculate agreement for Waking, NonREM, and REM sleep:
numberMismatchRed = zeros(6,4);
for i = 1:6
    numberMismatchRed(i,1:4) = [(numberMismatch(i,1) + numberMismatch(i,4) +numberMismatch(i,7)) (numberMismatch(i,2) + numberMismatch(i,6)) numberMismatch(i,3) numberMismatch(i,5)];
end

numberMismatchRed2 = [(numberMismatchRed(1,:) + numberMismatchRed(4,:)); (numberMismatchRed(2,:) + numberMismatchRed(6,:)); numberMismatchRed(3,:); numberMismatchRed(5,:)];
clear numberMismatchRed

totalAgreeRed = 0;

manualRedTotals = sum(numberMismatchRed2,2);


for i = 1:4
    manualVsAutoRed(i) = numberMismatchRed2(i,i)/manualRedTotals(i);
    totalAgreeRed = totalAgreeRed + numberMismatchRed2(i,i);
end

totalEpochRed = sum(manualRedTotals);
percentAgreeRed = totalAgreeRed/totalEpochRed;
%Write reduced agreement results:
xlswrite(resultExcelFileName, cellstr('Reduced % Agreement'), sheetName, 'A12');
xlswrite(resultExcelFileName, percentAgreeRed, sheetName, 'B12');

xlswrite(resultExcelFileName,cellstr('Waking'), sheetName, 'B14');
xlswrite(resultExcelFileName,cellstr('NonREM'), sheetName, 'C14');
xlswrite(resultExcelFileName,cellstr('REM'), sheetName, 'D14');
xlswrite(resultExcelFileName,cellstr('Unhooked'), sheetName, 'E14');
xlswrite(resultExcelFileName,cellstr('Manual Total'), sheetName, 'F14');
xlswrite(resultExcelFileName,cellstr('AS % Agree'), sheetName, 'G14');

xlswrite(resultExcelFileName,cellstr('Waking'), sheetName, 'A15');
xlswrite(resultExcelFileName,cellstr('NonREM'), sheetName, 'A16');
xlswrite(resultExcelFileName,cellstr('REM'), sheetName, 'A17');
xlswrite(resultExcelFileName,cellstr('UH'), sheetName, 'A18');

xlswrite(resultExcelFileName,[numberMismatchRed2 manualRedTotals manualVsAutoRed'],sheetName, 'B15');

%Calculate agreement for ALL sleep:
manualTotalSleep = manualRedTotals(2) + manualRedTotals(3);
sleepEpochMatch = numberMismatchRed2(2,2) + numberMismatchRed2(2,3) + numberMismatchRed2(3,2) + numberMismatchRed2(3,3);
asSleepAgreement = sleepEpochMatch/manualTotalSleep;

xlswrite(resultExcelFileName, cellstr('Sleep % Agreement'), sheetName, 'A20');
xlswrite(resultExcelFileName, asSleepAgreement, sheetName, 'B20');
end