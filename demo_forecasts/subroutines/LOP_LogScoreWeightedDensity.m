function [AggregatedDensity,WeightsMat] = LOP_LogScoreWeightedDensity(LogScoresMatrix,DensityMatrix,ForgettingFactor)

[TT,NoModels] = size(LogScoresMatrix);

% RECURSIVE LOG SCORE WEIGHTS 
% Cumulates log scores and at each point for all models examined. 

% Initialize at Equal Combination 
cl_scores = (1/NoModels)*ones(1,NoModels);

for tt=1:TT
    
    % Translate in Weights
    Weights=exp(cl_scores)./sum(exp(cl_scores));  
    WeightsMat(tt,:) = Weights; 
    
    density_t = DensityMatrix(:,:,tt);
    AggregatedDensity(:,tt)=density_t*Weights'; 
    
    l_scores = LogScoresMatrix(tt,:);
    
    %Update the weights: Cumulate (with possible discounting)
    w_s = ForgettingFactor; 
    cl_scores=(1-w_s)*cl_scores+w_s*l_scores; 
end
