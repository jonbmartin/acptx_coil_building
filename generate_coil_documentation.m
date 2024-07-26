% best solution
load test_grow_threshold/results_SAR_shift0_vflip0_hflip1.mat;
%load test/results_SAR_shift3_vflip0_hflip0.mat
load data/acptx_to_xin_array_convention_table.mat
% load data  
load data/head_2_b1field.mat
load data/mask_head_2.mat
b1 = Hx_allcoils + Hy_allcoils;
b1 = permute(b1,[2, 3, 4, 1]);
for coil_index = 1:30
    b1(:,:,:,coil_index) = b1(:,:,:,coil_index).* mask;
end
best_sol_index = 7;
slice_to_examine = 30;
mip_axis = 1;

%% actually construct plots

m_best = squeeze(results.all_m(best_sol_index,:,:,:));
b_best = squeeze(results.all_b(best_sol_index,:,:,:));

solution_coil_groupings = results.solution_coil_groupings;
solution_coil_groupings =  acptx_indexing_to_xin_indexing(solution_coil_groupings, xin_space_to_acptx_space);
solution_coil_index_mats = results.solution_coil_index_mats;
tiledlayout(8,5, "TileSpacing","compact")
net_b1 = zeros(61, 61, 61);
for grouping_index = 1:8 
    coil_group = solution_coil_groupings(grouping_index,:);
    group_b1 = zeros(size(net_b1));
    for coil_index = 1:4
        coil = coil_group(coil_index);
        if coil < 1
            % placeholder case - no coil
            nexttile, imagesc(zeros([61, 61]));
            xticks([]), yticks([])
            continue
        else
            coil_b1 = b1(:,:,:,coil)*b_best(coil,slice_to_examine);
            group_b1 = group_b1 + coil_b1;
            net_b1 = net_b1 + coil_b1;
            nexttile, imagesc(rot90(squeeze(max(abs(coil_b1),[],mip_axis)))), colormap jet, clim([0, 5]);
            xticks([]), yticks([])
        end

        % label plots
        if coil_index ==1
            ylabel(strcat('coil group ',int2str(grouping_index)))
        end
        title(strcat('|coil ',int2str(coil_index),' b1|'))

    end
    nexttile, imagesc(rot90(squeeze(max(abs(group_b1),[],mip_axis)))), colormap jet, clim([0, 5]);
    title('|GROUP b1|')
    xticks([]), yticks([])

end

% figure, imagesc(abs(net_b1)), title('|net b1|'), 
% colormap jet, clim([0, 1.5]),xticks([]), yticks([]);

