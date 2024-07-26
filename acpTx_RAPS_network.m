function network = acpTx_RAPS_network(weights,type,insertionLoss)

% Lays out a lossless 1-to-N array compression network that
% implements the attenuations in the weights vector,
% using ratio-adjustable power splitters (RAPS) (Yan, ISMRM 2017)

% type: layout, 'standardTree' (standard tree with 2^{level-1} splitters per
%       level, all outputs in the last level) or 'onePerStage'
%       (one splitter per stage, each stage has an output)

% If insertionLoss is provided (in dB), it accounts for that loss,
% since different coils will see Nstages * insertionLoss

% The output network array is Nstages x Nc,
% where rows correspond to stages, and columns to coils.

% It may be also possible to derive these weights using a polynomial
% analysis, where each output is one order.

% Written by Will Grissom, Vanderbilt University, 2017

weights = weights.^2; % the math is in terms of power

if nargin < 3
    insertionLossLinear = 1;
else
    insertionLossLinear = 10^(-insertionLoss/10); % divide by 10 since math is in terms of power
end

if nargin < 2
    type = 'standardTree';
end

switch type
    
    case 'onePerStage'
        
        Nc = length(weights);
        network = zeros(Nc-1,Nc);
        
        [~,sortInds] = sort(weights); % sort the weights, smallest-to-largest
        
        % Account for insertion loss, relative to first output. The n = 2nd
        % output through the n = (Nc-1)th output see n-1 attenuations, and the last
        % also sees Nc-2 attenuations.
        weights(flip(sortInds(2:end-1))) = ...
            weights(flip(sortInds(2:end-1)))./insertionLossLinear.^(1:Nc-2);
        weights(sortInds(1)) = weights(sortInds(1))/insertionLossLinear^(Nc-2);
        
        % start with the smallest two weights
        % solve equations: xi + xj = 1, xi/xj = weights(i)/weights(j)
        A = [1 1;1 -weights(sortInds(1))/weights(sortInds(2))];
        y = [1;0];
        x = A\y;
        network(end,sortInds(1)) = x(1);
        network(end,sortInds(2)) = x(2);
        
        for ii = 2:Nc-1 % loop over remaining stages, starting at the end
            A = [1 1;...
                network(end-(ii-2),sortInds(ii)) ...
                -weights(sortInds(ii))/weights(sortInds(ii+1))];
            x = A\y;
            network(end-(ii-2)-1,sortInds(ii)) = x(1);
            network(end-(ii-2)-1,sortInds(ii+1)) = x(2);
        end
        
    case 'standardTree'
        
        if insertionLossLinear ~= 1
            warning 'Insertion Loss compensation is ignored for Standard Tree, since it doesnt matter if number of outputs is even.'
        end
        
        if log2(length(weights)) - floor(log2(length(weights))) > 0 %rem(length(weights),2)
            % append zeros to get to closest power of 2
            fprintf(['Please note: since # coils is not a power of 2,\n' ...
                'we will append %d columns on the end of the network\n' ...
                'that can be disregarded.\n'],2^ceil(log2(length(weights)))-length(weights));
            weights = [weights(:).' eps*ones(1,2^ceil(log2(length(weights)))-length(weights))]; 
        end
        Nc = length(weights);
        
        % start at last (output) level, move backwards
        % for each consecutive output pair/splitter, solve equations  xi + xj = 1, xi/xj = weights(i)/weights(j)
        network = zeros(log2(Nc),Nc);
        y = [1;0]; % RHS
        for ii = 1:2:Nc % loop over splitters on this level
            A = [1 1;1 -weights(ii)/weights(ii+1)];
            x = A\y;
            network(end,ii) = x(1);
            network(end,ii+1) = x(2);
        end
        
        % loop over remaining stages
        for jj = log2(Nc)-1:-1:1 % loop over levels
            for ii = 1:2^(log2(Nc)-jj+1):Nc % loop over splitters on this level
                A = [1 1;...
                     prod(network(jj+1:end,ii))/prod(network(jj+1:end,ii+2^(log2(Nc)-jj)))...
                    -weights(ii)/weights(ii+2^(log2(Nc)-jj))];
                x = A\y;
                network(jj,ii) = x(1);
                network(jj,ii+2^(log2(Nc)-jj)) = x(2);
            end
        end
        
    otherwise
        
        error 'Unrecognized network topology. Valid options are standardTree and onePerStage.'
        
end

network = sqrt(network); % go back to voltages
