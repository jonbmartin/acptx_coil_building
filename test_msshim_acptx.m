% load the b1 maps
load data/LoopTx30_Head1_B1p_Axial.mat; % b1, mask
% load the coil-channel groups
load default_coil_groupings; % coil_groupings
coil_groupings(5, 3:4) = -1; % channel 5 has two repeated # 17's

[Nx, Ny, nSlices, nCoils] = size(b1);
nChannels = size(coil_groupings, 1);

% get the A matrices and initial phases for each slice 
Aall = {};
AinvAll = {};
phsAll = {};
lamda = 0.1; 
[x, y] = ndgrid(1 : Nx, 1 : Ny);
for ii = 1 : nSlices 

    tmp = permute(squeeze(b1(:, :, ii, :)), [3 1 2]);
    tmp = tmp(:, :).';
    tmp = tmp(logical(mask(:, :, ii)), :);
    Aall{ii} = tmp; 
    AinvAll{ii} = inv(tmp'*tmp + lamda * eye(nCoils)) * tmp'; 

    % find the slice's centroid
    xMean = round(mean(x(col(logical(mask(:, :, ii))))));
    yMean = round(mean(y(col(logical(mask(:, :, ii))))));

    b(:, ii) = squeeze(exp(-1i*angle(b1(xMean, yMean, ii, :))));

    phsAll{ii} = angle(tmp * b(:, ii));


end

% for debug - get all the shim patterns at this point 
minit = [];
for ss = 1 : nSlices 
    minit(:, :, ss) = embed(Aall{ss} * b(:, ss), logical(mask(:, :, ss)));
end

% % initial shim - this gives a nice result but does not help the later svd-thresholded run, so lets skip
% nIters = 100;
% cost = zeros(nIters, 1);
% b = ones(nCoils, nSlices);
% for ii = 1 : nIters 

%     % update the shims 
%     for ss = 1 : nSlices 
%         b(:, ss) = AinvAll{ss} * exp(1i * phsAll{ss});
%     end 

%     % update phases
%     for ss = 1 : nSlices 
%         phsAll{ss} = angle(Aall{ss} * b(:, ss));
%     end 

%     % measure cost 
%     for ss = 1 : nSlices 
%         cost(ii) = cost(ii) + sum(abs(exp(1i * phsAll{ss}) - Aall{ss} * b(:, ss)).^2);
%     end
%     cost(ii) = cost(ii) + lamda * sum(abs(b(:)).^2);

%     printf('Iteration %d. Cost: %0.1f', ii, cost(ii));

% end

% % for debug - get all the shim patterns at this point 
% m = [];
% for ss = 1 : nSlices 
%     m(:, :, ss) = embed(Aall{ss} * b(:, ss), logical(mask(:, :, ss)));
% end 

nIters = 500; 
nCgIters = 3;
cost = zeros(nIters, 1);
b = ones(nCoils, nSlices);
for ii = 1 : nIters 

    % update the shims 
    for ss = 1 : nSlices 
        [xS, info] = qpwls_pcg(b(:,ss),[Aall{ss};sqrt(lamda)*eye(size(b(:,ss),1))],1,...
                                [exp(1i*phsAll{ss});0*b(:,ss)],0,0,1,nCgIters,ones(size(b,1),1));
        b(:, ss) = xS(:, end);
        %b(:, ss) = AinvAll{ss} * exp(1i * phsAll{ss});
    end 

    % for debug - get all the shim patterns at this point 
    m = [];
    for ss = 1 : nSlices 
        m(:, :, ss) = embed(Aall{ss} * b(:, ss), logical(mask(:, :, ss)));
    end
    figure(1);im(121,m);title 'before trunc',colorbar

    % truncate SVDs
    sAll = zeros(4, 1); % to keep track of sizes of singular values
    for cc = 1 : nChannels
        coils_chan = coil_groupings(cc, :);
        coils_chan = coils_chan(coils_chan ~= -1);
        % pull the shims for this channel into a matrix
        bSub = b(coils_chan, :);
        % svd of the matrix
        [u, s, v] = svd(bSub, 'econ');
        % truncate or soft-threshold all the singular values except at (1,1)
        s(2:end) = max(0, s(2:end) - 0.25);
        %st = st * sum(s(:)) / sum(st(:));
        %s = st;
        %s(2 : end) = 0;
        % reconstitute the matrix
        bSub = u * s * v';

        % put it back with the rest
        b(coils_chan, :) = bSub;
        sAll(1:length(coils_chan)) = sAll(1:length(coils_chan)) + diag(s);
    end

    % [u, s, v] = svd(b, 'econ');
    % st = zeros(size(s));
    % st(1:8, 1:8) = s(1:8, 1:8);
    % b = u * st * v';

    % update phases
    for ss = 1 : nSlices 
        phsAll{ss} = angle(Aall{ss} * b(:, ss));
    end 

    % for debug - get all the shim patterns at this point 
    m = [];
    for ss = 1 : nSlices 
        m(:, :, ss) = embed(Aall{ss} * b(:, ss), logical(mask(:, :, ss)));
    end
    figure(1);im(122,m);title 'after trunc', colorbar, sAll, drawnow %, pause; 



    % measure cost 
    for ss = 1 : nSlices 
        cost(ii) = cost(ii) + sum(abs(exp(1i * phsAll{ss}) - Aall{ss} * b(:, ss)).^2);
    end
    cost(ii) = cost(ii) + lamda * sum(abs(b(:)).^2);

    printf('Iteration %d. Cost: %0.1f', ii, cost(ii));

end

