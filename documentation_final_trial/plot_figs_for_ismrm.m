% figure 5
close all
tiledlayout(1,2,'Padding', 'compact', 'TileSpacing', 'compact')

% show CP mode shim 
load("../data/head_1_b1field.mat")
%H = sqrt(Hx_allcoils.^2 + Hy_allcoils.^2);
H = Hx_allcoils + 1j * Hy_allcoils;

phs = angle(H(:,30,30,30));
phs_wt = exp(-phs*1j);

n_slices = 61;
cp_shimmed_head = zeros(61,61,61);
load("../data/mask_head_1.mat")
for ii=1:n_slices
    phs = angle(H(:,30,30,ii));
    phs_wt = exp(-phs*1j);
    slice_shim = squeeze(sum(phs_wt.*H(:,:,:,ii),1));
    cp_shimmed_head(:,:,ii) = squeeze(mask(:,:,ii)).*slice_shim;
end

load("results_SAR_shift6_vflip0_hflip1.mat")

ex_1 =results.all_m(9,:,:,1:4:48);
ex_2 =results.all_m(6,:,:,1:4:48);

for ii = 1:size(ex_1,4)
    ex_1(:,:,:,ii) = rot90(squeeze(ex_1(:,:,:,ii)),-1);
    ex_2(:,:,:,ii) = rot90(squeeze(ex_2(:,:,:,ii)),-1);
end
cov_ex_1 = std(nonzeros(abs(ex_1)))/mean(nonzeros(abs(ex_1)));
cov_ex_2 = std(nonzeros(abs(ex_2)))/mean(nonzeros(abs(ex_2)));

nexttile,
im(ex_1)
text(-0.15,0.95,'a)','Units','normalized','FontSize',12)
title('')
xticklabels(''), yticklabels('')
clim([0, 1.3])
colormap turbo
h = colorbar;
h.Label.String = "|B_{1}^{+}| (a.u.)";
h.Label.Rotation = 270;
h.Label.VerticalAlignment = "bottom";

%load("results_SAR_shift6_vflip0_hflip1.mat")
nexttile, 
im(ex_2)
text(-0.15,0.95,'b)','Units','normalized','FontSize',12)
title('')
xticklabels(''), yticklabels('')
clim([0, 1.3])

colormap turbo
h = colorbar;
h.Label.String = "|B_{1}^{+}| (a.u.)";
h.Label.Rotation = 270;
h.Label.VerticalAlignment = "bottom";
