% % load the b1 maps, set up directory info
% load data/LoopTx30_Head1_B1p_Axial.mat; % b1, mask
% load data/LoopTx8_Head1_Cond3dxyz.mat; % cond
% load data/LoopTx8_Head1_MassDensity3dxyz.mat; % mass density
% load data/LoopTx30_Head1_ExEyEz3dxyz.mat; % individual e fields
% load data/acptx_to_xin_array_convention_table.mat
% 
% % format into structure
% maps.Ex_out = Ex_z_x_y_all;
% maps.Ey_out = Ey_z_x_y_all;
% maps.Ez_out = Ez_z_x_y_all;
% Emag = sqrt(Ex_z_x_y_all.^2 + Ey_z_x_y_all.^2 + Ez_z_x_y_all.^2);
% maps.cond3D_out = Cond3dxyz;
% maps.dens3D_out = MassDensity3dxyz;
% 
% S = calc_10g_sar(maps);
% Sglobal = sum(S,3)/size(S,3);
% R = full(buildSARReg(Sglobal,1)); % build SAR regularization mat
% % get to right dimension
% B1 = b1;
% datasize = 101;
% n_slices_total = size(B1,3);


load data/acptx_to_xin_array_convention_table.mat

load data/head_1_b1field.mat
load data/mask_head_1.mat
load data/head_1_global_sar.mat
B1 = (Hx_allcoils + Hy_allcoils);
mask_coilcopy = permute(repmat(mask,1,1,1,30),[4,1,2,3]);
B1 = B1.* mask_coilcopy;
B1 = permute(B1,[2,3,4,1]); % x, y, slice, coil
bottom_slice = 10;
max_nonzero_slice = 54;
datasize = 61;
B1_1 = B1(:,:,bottom_slice:max_nonzero_slice,:);
mask_1 = mask(:,:,bottom_slice:max_nonzero_slice);
Sglobal_1 = S;

load data/head_2_b1field.mat
load data/mask_head_2.mat
load data/head_2_global_sar.mat
B1 = (Hx_allcoils + Hy_allcoils);
mask_coilcopy = permute(repmat(mask,1,1,1,30),[4,1,2,3]);
B1 = B1.* mask_coilcopy;
B1 = permute(B1,[2,3,4,1]); % x, y, slice, coil
bottom_slice = 10;
max_nonzero_slice = 54;
datasize = 61;
B1_2 = B1(:,:,bottom_slice:max_nonzero_slice,:);
mask_2 = mask(:,:,bottom_slice:max_nonzero_slice);
Sglobal_2 = S;

%Concatenate
concatenate = true;
if concatenate
    Sglobal= Sglobal_1 + Sglobal_2;
    B1 = cat(3,B1_1,B1_2);
    mask = cat(3, mask_1, mask_2);
    n_slices_total = size(B1,3);
else
    Sglobal = Sglobal_2;
    B1 = B1_2;
    mask = mask_2;
    n_slices_total = size(B1,3);
end

%%
local_savedir = 'test_grow_threshold';

% load the coil-channel groups
load default_coil_groupings; % coil_groupings
coil_groupings(5, 3:4) = -1; % channel 5 has two repeated # 17's

% build the indexing matrix - matches the diagram
coil_index_mat = build_array(30, 3, 10);

% apply the operation to the matrix
n_shifts = 9; vertical_flip = [false, true]; horizontal_flip = [false, true];
total_num_configs = (n_shifts+1) * 2 * 2;

% regularization v[alues to test
n_reg_vals = 3; % formerly 15
lamda_v = logspace(-1, 1, n_reg_vals);% works ok but overregularizes SAR

lamda_v = logspace(-3, -1, n_reg_vals);
%lamda_v = logspace(-3, -1, n_reg_vals);
%lambda_v = ones(size(lamda_v))*0.025;
% general optimization parameterst
nIters = 500; % 500
nCgIters = 3;
plotting = false;


% TODO: want to create some cell storate array to store L-curve, coil
% config for each coil grouping combination 

best_config = [];
best_coil_index_mat = [];
best_b = [];

lcurve_results = zeros(3, n_reg_vals);%1 = total, 2 = sar, 3 = rmse

% set up structure to store results - is overwritten with each config
% tested and then saved to a unique mfile
results.lamda_v = lamda_v;
results.nIters = nIters;
results.nCgIters = nCgIters;

results.solution_coil_index_mats = zeros(3,10);
results.solution_coil_groupings = zeros(8,4);
results.shift = [];
results.hflip = [];
results.vflip = [];
results.solution_lcurve = [];
results.all_m = zeros(n_reg_vals, datasize, datasize, n_slices_total);
results.all_b = zeros(n_reg_vals, 30, n_slices_total);

% loop through all candidate configs, and perform an optimization across
% regularization values.
for shift_index = 0: 9
    for hflip_index = 1:length(horizontal_flip)
        for vflip_index = 1:length(vertical_flip)
            results.shift = shift_index;
            results.hflip = horizontal_flip(hflip_index);
            results.vflip = vertical_flip(vflip_index);

            % create a new config
            coil_index_mat = build_array(30, 3, 10);
            [coil_index_mat, new_coil_grouping] = output_new_config(coil_index_mat, shift_index, vertical_flip(vflip_index), horizontal_flip(hflip_index));
            new_coil_grouping_xin_convention = acptx_indexing_to_xin_indexing(new_coil_grouping, flipud(xin_space_to_acptx_space));
            results.solution_coil_groupings = new_coil_grouping;
            results.solution_coil_index_mats = coil_index_mat;
            results.solution_coil_groupings_xin_space = new_coil_grouping_xin_convention;
            
            % make sure that config is indexing the proper b1 coils
            

            for lamda_index = 1:n_reg_vals
                lamda = lamda_v(lamda_index);
                
                %pretty print the iteration's info
                disp(['shift =  ', num2str(shift_index)]);
                disp(['horizontal flip =  ', num2str(horizontal_flip(hflip_index))]);
                disp(['vertical flip =  ', num2str(vertical_flip(vflip_index))]);
                disp(['lamda value # ',num2str(lamda_index),', = ', num2str(lamda)]);

                % perform the optimization
                [cost, b, m] = msshim_acptx_fx(B1, mask, lamda, Sglobal,...
                    new_coil_grouping_xin_convention, nIters, nCgIters, plotting);
                
                % store the cost
                [min_total_cost, min_index] = min(cost{1,3});
                min_index = length(cost{1,3});
                lcurve_results(1,lamda_index) = min_total_cost;
                sar_cost = cost{1,2};
                lcurve_results(2,lamda_index) = sar_cost(min_index);
                rmse_cost = cost{1,1};
                lcurve_results(3,lamda_index) = rmse_cost(min_index);

                % store the shim and m
                results.all_m(lamda_index, : ,: ,:) = m;
                results.all_b(lamda_index, :, :) = b;
            end
            % store l curve results for the given configuration
            results.rmse_cost_v = lcurve_results(3,:);
            results.sar_cost_v = lcurve_results(2, :);
            
            config_string = strcat('shift',num2str(shift_index), ...
                '_vflip',num2str(results.vflip),'_hflip',num2str(results.hflip))
            disp('==========================================')
            disp(strcat('Config ',config_string, ' completed, saving ...'))
            disp('==========================================')
        
            % save every configs results - messy but necessary b/c of file size
            savestring = strcat(local_savedir, '/results_SAR_',...
                config_string, '.mat');
            save(savestring, 'results')
        end
    end

end


