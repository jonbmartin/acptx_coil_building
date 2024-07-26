D = '/home/jonathan/PycharmProjects/minibatch_shimming/acptx/acptx_coil_building/test_grow_threshold/';

S = dir(fullfile(D,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.


for ii = 3:numel(S) % 3 to skip '.' and '..'
    file = strcat(D,S(ii).name);
    disp(file)
    load(file);

    m = results.all_m;
    [nreg,datasizex, datasizey, nslices] = size(m);
    hole_mask = zeros(nreg,1);

    for reg_index = 1:nreg


        
        mslice = squeeze(m(reg_index,:,:,:));
        mslice = abs(mslice);
        mask = mslice > 0;
        voids = mslice < 0.2;  
        null = sum(sum(sum(mask .* voids)));
        if null > 0
            hole_mask(reg_index,:) = 1;
        end
        
    end
    sum(hole_mask)
end