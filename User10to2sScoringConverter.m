function [manualScored2sEpochs,scoredTS, fileLabel] = User10to2sScoringConverter
%Request user input to select manually scored file:
working_dir = pwd;
current_dir = 'C:\SleepData\Results'; % This is the default directory that opens.
cd(current_dir);
[filename, pathname] = uigetfile('*.xls', 'Select manually scored file'); %This waits for user input and limits selection to .xls files.

% Check for whether or not a file was selected
if isempty(filename) || isempty(pathname)
    uiwait(errordlg('You need to select a manually scored file. Please try again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    userFile= fullfile(pathname, filename);
end
fileLabel = strrep(filename, '.xls', '');
% Read in the manually scored file (10s epochs) in Excel .XLS format:
manualScoredEpochs = xlsread(userFile);
calcEpochSize = manualScoredEpochs(2,2) - manualScoredEpochs(1,2);
if calcEpochSize > 9
    num10sEpochs = size(manualScoredEpochs,1); %Find # of 10s epochs in manually scored file.
    manualScored2sEpochs = zeros(num10sEpochs*5,1); %Create a vector variable with appropriate length for 2s epoch conversion.
    scoredTS = manualScoredEpochs(:,2);
    %Create five 2s epochs for each 10s epoch with the same stage:
    for i=1:num10sEpochs
        manualScored2sEpochs((5*i-4):(5*i)) = manualScoredEpochs(i,3);
    end
elseif calcEpochSize < 3
    manualScored2sEpochs(:,1) = manualScoredEpochs(:,3);
else
    uiwait(errordlg('The selected manually scored file is not in 10s or 2s epochs. Please try again.',...
        'ERROR','modal'));
    %Need to add user choice dialog to try again or terminate the program
end
