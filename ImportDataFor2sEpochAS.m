function[PSDValues]=ImportDataFor2sEpochAS(Thresholds,filename1,filename2,filename3,varargin)

if length(varargin) == 3
    filename4=char(varargin(1,1));
    lowertimestamp = str2num(char(varargin(1,end-1)));
    uppertimestamp = str2num(char(varargin(1,end)));
%     lowerbound=str2num(char(varargin(1,end-1)));
%     upperbound=str2num(char(varargin(1,end)));
else
    lowertimestamp = str2num(char(varargin(1,end-1)));
    uppertimestamp = str2num(char(varargin(1,end)));
%     lowerbound=str2num(char(varargin(1,end-1)));
%     upperbound=str2num(char(varargin(1,end)));
    filename4=[];
end

cwd = pwd;
cd(tempdir);

cd(cwd);

nsamp = []; % Number of samples per bin of time


% Neuralynx System
waitbar(0.1,waithandle,'Converting EMG from Neuralynx CSC to Matlab data ...');
figure(waithandle),pause(0.2),
[Timestamps,SF,Samples] = Nlx2MatCSC(filename1,[1 0 1 0 1],0,4,[lowertimestamp uppertimestamp]);
samples = double(Samples(:)');
clear Samples
SampFreq1 = SF(1);
if boundIndex == 1
    DS = (1:1:10);
    DSampSF = SampFreq1./DS;
    indSampfactor = find(DSampSF >= 250);
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
waitbar(0.4,waithandle,'Filtering the EMG data ...'); 
figure(waithandle),pause(0.2),
physInput = 1;  %Needed to select proper error box in HeaderADBit.
ADBit2uV = HeaderADBit(filename1, physInput);    %Calls a function to extract the AD Bit Value.
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
 
clear filtered_samples samples adfreq

%  ******    EEG FILE extraction   *********

% Neuralynx System
    waitbar(0.6,waithandle,' Converting EEG from Neuralynx CSC to Matlab data ...');
    figure(waithandle),pause(0.2)
    [Timestamps,SF,Samples]=Nlx2MatCSC(filename2,[1 0 1 0 1],0,4,[lowertimestamp uppertimestamp]);
    samples=double(Samples(:)');
    clear Samples
    SampFreq2=SF(1);
    Fs=round(SampFreq2/sampfactor);
    waitbar(0.8,waithandle,'Filtering the EEG data ...'); 
    figure(waithandle),pause(0.2),
    [EEG_TIMESTAMPS,EEG_SAMPLES] = generate_timestamps_from_Ncsfiles(Timestamps,samples, exactLow, exactHi, nsamp);
    physInput = 2;  %Needed to select proper error box in HeaderADBit.
    ADBit2uV = HeaderADBit(filename2, physInput);    %Calls a function to extract the AD Bit Value.
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

EPOCHtime=[]; EPOCHnum=[];
if isempty(filename4)==0    %Checks to see if Training file exists
    try
        [t_stamps]=xlsread(filename4);  %Imports Training file data
    catch
        uiwait(errordlg('Check if the file is saved in Microsoft Excel format.',...
            'ERROR','modal'));
    end
    EPOCHtime = t_stamps(1:end,2);  % Makes column of downsampled time steps new vector 
    EPOCHnum = t_stamps(1:end,3);   % Makes column of scored states for each time step above
    clear t_stamps
end
fprintf('          Calculating the FFT and getting the power values....\n');
% This is to find the EPOCHSIZE of 10 sec
index=find((EMG_TIMESTAMPS(1)+9.999 < EMG_TIMESTAMPS) & (EMG_TIMESTAMPS < EMG_TIMESTAMPS(1)+10.001));
if(isempty(index)) == 1
    index=find((EMG_TIMESTAMPS(1)+9.99 < EMG_TIMESTAMPS) & (EMG_TIMESTAMPS < EMG_TIMESTAMPS(1)+10.01));
end
diff= EMG_TIMESTAMPS(index(1):index(end)) - (EMG_TIMESTAMPS(1)+10);
[minimum,ind]=min(abs(diff));
try
    EPOCHSIZE=index(ind);
catch
        fprintf('There is an error in calculating the EPOCHSIZE of 10sec in read_n_extract_datafiles\n');
end
%Take the fft of the entire timeframe we have, calcuate power in bins
[PSDValues]=fft_psd_and_statescore_of_epoch(Thresholds);
waitbar(1,waithandle, 'Finished converting.. Now Loading the data ..');
figure(waithandle), pause(0.3);
close(waithandle);