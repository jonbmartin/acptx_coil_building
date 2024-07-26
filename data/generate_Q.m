load head_1_efield.mat
load sarfact_head_1.mat

% sarfact = sarfact_head_1;
Nc = 30;
for coil_index = 1:Nc
    Ex_allcoils(coil_index,:,:,:) = squeeze(Ex_allcoils(coil_index,:,:,:)).*abs(sarfact); 
    Ey_allcoils(coil_index,:,:,:) = squeeze(Ey_allcoils(coil_index,:,:,:)).*abs(sarfact); 
    Ez_allcoils(coil_index,:,:,:) = squeeze(Ez_allcoils(coil_index,:,:,:)).*abs(sarfact); 
end

% flatten the spatial dimensions
Ex_allcoils = Ex_allcoils(:,:,:);Ex_allcoils = Ex_allcoils(:,:);
Ey_allcoils = Ey_allcoils(:,:,:);Ey_allcoils = Ey_allcoils(:,:);
Ez_allcoils = Ez_allcoils(:,:,:);Ez_allcoils = Ez_allcoils(:,:);

Ns = size(Ex_allcoils,2);

S = zeros(Nc,Nc,Ns);
for spatial_index = 1:Ns
    S(:,:,spatial_index) = Ex_allcoils(:,spatial_index)*Ex_allcoils(:,spatial_index)' + Ey_allcoils(:,spatial_index)*Ey_allcoils(:,spatial_index)' + Ez_allcoils(:,spatial_index)*Ez_allcoils(:,spatial_index)';
end
S = mean(S,3);