% State vectors:
tenSecStates = [1 1 1 0 1]; %This is setting up the history for start of 
tenSecSleep = [0 0 0 0 0]; %the program for the 10 seconds prior to programs start.
%real3EpochHistory = [0 0 0];  % This information is needed in the auto-scoring algorithm.
autoScoredEpochs = [];  % Record the auto-scored epochs and save at end of program.

for i= 1:numOfEpochs
    % FULL AUTOSCORER LOGIC
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
        if sigmaPower(i) > (Thresholds(6) + 2*Thresholds(7))  % This line is in the analysis of 1 s epochs.
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
        if isequal(tenSecStates(3:5), [1 1 1])
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

%This is the END of the auto-scoring component.  Now set the sleep/wake
%value in the 'tenSecStates' vector for the air puff logic.
    % Determine state via logic.
    % State = 0 is asleep, State = 1 is awake
    % Push-Pop states.
    switch tempState
        case 1
            tenSecStates = [tenSecStates(2:5) 1]; %Scored as AWAKE.
        case 2
            tenSecStates = [tenSecStates(2:5) 0]; %Scored as ASLEEP.
        case 3
            tenSecStates = [tenSecStates(2:5) 0]; %Scored as ASLEEP.
        case 4
            tenSecStates = [tenSecStates(2:5) 1]; %Scored as AWAKE.
        case 5
            tenSecStates = [tenSecStates(2:5) 1]; %Scored as AWAKE.
        case 6
            tenSecStates = [tenSecStates(2:5) 0]; %Scored as ASLEEP.
        % case 7 will not occur since the epoch will always be
        % scored.
        case 8
            tenSecStates = [tenSecStates(2:5) 1]; %Scored as AWAKE.
    end
end