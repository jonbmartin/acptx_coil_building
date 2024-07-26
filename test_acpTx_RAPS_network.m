% Script to test the acpTx_RAPS_network function
% Written by Will Grissom, Vanderbilt University, 2017

%% Define the coil weights to be implemented by the network
% weights must all be positive!
weights = [0.13 0.55 0.2 1]; % Channel 2 network from Yan, MRM 2016, Fig. 6
%weights = rand(1,4); % random between 0 and 1


%% Compute both network topologies that implement the weights
% The computed 'network' matrices are Nstages x Nc,
% where rows correspond to stages, and columns to coil outputs.
% For a standardTree network, non-zero elements in a row that are above
% non-zero elements in the next row correspond to the connections to that
% stage. All outputs are in the last stage.
% For a onePerStage network, non-zero elements in a row that are above
% non-zero elements in the next row correspond to the connections to that
% stage. Non-zero elements in a row that are above zero elements in the 
% next row correspond to the output to that column's coil. 
%insertionLoss = 0.3; % insertion loss in dB
insertionLoss = 0; % insertion loss in dB
type = 'standardTree';
networkStandardTree = acpTx_RAPS_network(weights,type,insertionLoss);

type = 'onePerStage';
networkOnePerStage = acpTx_RAPS_network(weights,type,insertionLoss);

%% Plot the networks
subplot(121)
[Nstages,Ncoils] = size(networkStandardTree);
standardTreeNodes = [0 kron(1:length(weights)-1,[1 1])];
treeplot(standardTreeNodes);
axis square
[x,y] = treelayout(standardTreeNodes);
text(x(1),y(1),'  Input','FontSize',12);
netWeightsCompressed = networkStandardTree.';
netWeightsCompressed = netWeightsCompressed(:);
netWeightsCompressed = netWeightsCompressed(netWeightsCompressed > 0);
for ii = 1:length(netWeightsCompressed)-Ncoils
    text(x(ii+1),y(ii+1),['  ' num2str(netWeightsCompressed(ii))],...
        'FontSize',12);
end
for ii = length(netWeightsCompressed)-Ncoils+1:length(netWeightsCompressed)
    text(x(ii+1),y(ii+1),['  Coil ' num2str(ii-(length(netWeightsCompressed)-Ncoils)) ...
        ':\newline  ' num2str(netWeightsCompressed(ii))],...
        'FontSize',12);
end
axis off
title('Standard Tree','FontSize',16)

subplot(122)
[Nstages,Ncoils] = size(networkOnePerStage);
onePerStageNodes = zeros(1,2^(length(weights)-1)-1);
prevParentInd = 1; % the initial parent
parentInd = 2;
outputInd = 3; % the initial output
weightLabels = {};
weightLabels{1} = '  Input';
for ii = 1:Nstages-1
    % find the output from this stage
    outputWtInd = find(networkOnePerStage(ii,:) > 0 & networkOnePerStage(ii+1,:) == 0);
    %outputIndNew = parentInd+2;
    onePerStageNodes(outputInd) = prevParentInd;
    weightLabels{outputInd} = ...
        ['  Coil ' num2str(outputWtInd) ':\newline  ' num2str(networkOnePerStage(ii,outputWtInd))];
    % find the input to the next stage; it connects to previous parent and
    % is parent to next stage
    parentWtInd = find(networkOnePerStage(ii,:) > 0 & networkOnePerStage(ii+1,:) > 0);
    %parentIndNew = parentInd+1;
    weightLabels{parentInd} = ['  ' num2str(networkOnePerStage(ii,parentWtInd))];
    onePerStageNodes(parentInd) = prevParentInd;
    prevParentInd = parentInd;
    parentInd = outputInd+1;
    outputInd = outputInd+2;
end
% the two outputs from the last stage come from the last parent
outputInd = find(networkOnePerStage(end,:) > 0);
onePerStageNodes(end-1) = prevParentInd;
onePerStageNodes(end) = prevParentInd;
weightLabels{end+1} = ...
        ['  Coil ' num2str(outputInd(1)) ':\newline  ' num2str(networkOnePerStage(end,outputInd(1)))];
weightLabels{end+1} = ...
        ['  Coil ' num2str(outputInd(2)) ':\newline  ' num2str(networkOnePerStage(end,outputInd(2)))];
treeplot(onePerStageNodes);
axis square
[x,y] = treelayout(onePerStageNodes);
for ii = 1:length(x)
    text(x(ii),y(ii),weightLabels{ii},'FontSize',12);
end
axis off
title('One-Per-Stage','FontSize',16)

%% Validate the networks by calculating the overall weight for each coil
% from the network matrices

% Standard Tree
% loop through stages, multiplying each level by weights above it
[Nstages,Ncoils] = size(networkStandardTree);
networkIter = networkStandardTree;
for ii = 1:Nstages-1
    % duplicate weights into non-zero entries
    for jj = 1:Ncoils/(2^ii):Ncoils
        networkIter(ii,jj+1:jj+Ncoils/(2^ii)-1) = networkIter(ii,jj);
    end
end
weightsConstructed = prod(networkIter,1);
% truncate constructed weights if not power of 2
weightsConstructed = weightsConstructed(1:length(weights));

% compare to target weights after normalization
disp 'Standard Tree Test Error (%):'
100*norm(weights./max(weights) - weightsConstructed./max(weightsConstructed))/norm(weights./max(weights))

% One Per Stage
Ncoils = length(weights);
Nstages = Ncoils-1;
weightsConstructed = zeros(1,Ncoils);
networkIter = networkOnePerStage;
insertionLossLinear = 10^(-insertionLoss/20);
for ii = 1:Nstages-1
    
    % find the entry that terminates here; put it in that coil's weight
    weightsConstructed(networkIter(ii,:) > 0 & networkIter(ii+1,:) == 0) = ...
        networkIter(ii,networkIter(ii,:) > 0 & networkIter(ii+1,:) == 0);
    
    % find the entry that feeds the next stage, and multiply it downstream
    % into the rest of the network
    nextStageConnection = find(networkIter(ii,:) > 0 & ...
        networkIter(ii+1,:) > 0);
    networkIter(ii+1:end,:) = networkIter(ii+1:end,:)*...
        networkOnePerStage(ii,nextStageConnection)*insertionLossLinear;
    
end
weightsConstructed(networkIter(end,:) > 0) = networkIter(end,networkIter(end,:) > 0);

% compare to target weights after normalization
disp 'One-Per-Stage Test Error (%):'
100*norm(weights./max(weights) - weightsConstructed./max(weightsConstructed))/norm(weights./max(weights))
