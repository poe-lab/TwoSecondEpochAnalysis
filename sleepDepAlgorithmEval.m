function sleepDepAlgorithmEval(deltaPowerThreshold)

numOfEpochs = length(scored2sEpochs);
matrixComparison = zeros(2,8);
for i=1:numOfEpochs
    if deltaPower(i) > deltaPowerThreshold
%         autoState = 0;
        switch scored2sEpochs(i)
            case 1 %'AW'
                matrixComparison(1,1) = matrixComparison(1,1) +1;
            case 2 %'QS'
                matrixComparison(1,2) = matrixComparison(1,2) +1;
            case 3 %'RE'
                matrixComparison(1,3) = matrixComparison(1,3) +1;
            case 4 %'QW'
                matrixComparison(1,4) = matrixComparison(1,4) +1; 
            case 5 %'UH'
                matrixComparison(1,5) = matrixComparison(1,5) +1;
            case 6 %'TR'
                matrixComparison(1,6) = matrixComparison(1,6) +1;
            case 7 %'NS'
                matrixComparison(1,7) = matrixComparison(1,7) +1;
            case 8 %'IW'
                matrixComparison(1,8) = matrixComparison(1,8) +1;
        end 
    else
%         autoState = 1;
        switch scored2sEpochs(i)
            case 1 %'AW'
                matrixComparison(2,1) = matrixComparison(2,1) +1;
            case 2 %'QS'
                matrixComparison(2,2) = matrixComparison(2,2) +1;
            case 3 %'RE'
                matrixComparison(2,3) = matrixComparison(2,3) +1;
            case 4 %'QW'
                matrixComparison(2,4) = matrixComparison(2,4) +1; 
            case 5 %'UH'
                matrixComparison(2,5) = matrixComparison(2,5) +1;
            case 6 %'TR'
                matrixComparison(2,6) = matrixComparison(2,6) +1;
            case 7 %'NS'
                matrixComparison(2,7) = matrixComparison(2,7) +1;
            case 8 %'IW'
                matrixComparison(2,8) = matrixComparison(2,8) +1;
        end 
    end
end

totalEachState = sum(matrixComparison);

totalSleepEpochs = totalEachState(2) + totalEachState(3) + totalEachState(6);
totalSleepScoredAsSleep = matrixComparison(1,2) + matrixComparison(1,3) + matrixComparison(1,6);
percentSleepCorrect = 100*totalSleepScoredAsSleep/totalSleepEpochs;

totalWakeEpochs = totalEachState(1) + totalEachState(4) + totalEachState(8);
totalWakeScoredAsWake = matrixComparison(2,1) + matrixComparison(2,4) + matrixComparison(2,8);
percentWakeCorrect = 100*totalWakeScoredAsWake/totalWakeEpochs;

overallPerformance = 100* (totalSleepScoredAsSleep + totalWakeScoredAsWake)/(totalSleepEpochs + totalWakeEpochs);
percentOfStateAsSleep = 100 * matrixComparison(1,:)./totalEachState;

