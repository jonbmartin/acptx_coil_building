D = '/home/jonathan/PycharmProjects/minibatch_shimming/acptx/acptx_coil_building/test_grow_threshold/';

S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.

min_err_sol = [];
min_err = inf;
figure,
for ii = 3:numel(S) % 3 to skip '.' and '..'
    file = strcat(D,S(ii).name);
    disp(file)
    load(file);
    hold on, plot(results.rmse_cost_v,results.sar_cost_v)
    hold on, scatter(results.rmse_cost_v,results.sar_cost_v)
    min_total_config_err = min(results.sar_cost_v + results.rmse_cost_v);
    min_rmse_err = min(results.rmse_cost_v);
    disp(strcat('solution error = ',num2str(min_total_config_err)))
    disp(strcat('rmse error = ',num2str(min_rmse_err)))

    if min_total_config_err < min_err
        min_err = min_total_config_err;
        disp(strcat('min overall solution = ', S(ii).name, '; min_err = ', num2str(min_err)))
    end

end

xlabel('Error cost')
ylabel('SAR cost')