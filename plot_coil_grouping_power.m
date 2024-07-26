%load results_power/results_shift2_vflip0_hflip0.mat
load test_grow_threshold/results_SAR_shift0_vflip0_hflip1.mat;
load data/acptx_to_xin_array_convention_table.mat;
reg_index = 7;

best_sol_mag = abs(squeeze(results.all_b(reg_index,:,:)));

%do NOT want this in xin space, it's in acptx space
coil_groupings_mat = acptx_indexing_to_xin_indexing(results.solution_coil_groupings, xin_space_to_acptx_space);
%coil_groupings_mat = results.solution_coil_groupings;
normalize_weights = true;
coil_group_weight_mat = zeros(8,4);
for grouping_index = 1:8
    coil_group = coil_groupings_mat(grouping_index,:);
    coil_group = coil_group(coil_group>0); % prune negative placeholders
    coil_group_b = best_sol_mag(coil_group,:);
    coil_group_weights = max(coil_group_b,[],2);
    if normalize_weights
        coil_group_weights = coil_group_weights / max(coil_group_weights); % normalize
    end
    if size(coil_group_weights,1) == 2
        coil_group_weights = [coil_group_weights; 0; 0];
    end
    coil_group_weight_mat(grouping_index, :) = coil_group_weights;
end

% plot the relative powers for ALL slices together 
best_b = squeeze(results.all_b(reg_index, :, :));
best_b_ordered = [];
for grouping_index = 1:8
    coil_group = coil_groupings_mat(grouping_index,:);
    coil_group = coil_group(coil_group>0);
    best_b_ordered = [best_b_ordered; best_b(coil_group,:)];
end

figure,
imagesc(coil_group_weight_mat), colormap jet, xlabel('coil index'), ylabel('coil group'),
colormap
