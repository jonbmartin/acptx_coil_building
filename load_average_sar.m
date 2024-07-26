pwd
% do E-field
datadirs = {'head_1_efield', 'head_2_efield'};
for dataset_index = 1:length(datadirs)
    % load all of the data in the directory
    myDir = strcat(pwd,'/data/',datadirs{dataset_index}); %gets directory
    myFiles = dir(fullfile(myDir,'*')); %gets all wav files in struct
    
    for data_index = 1:size(myFiles,1)
        data_index
        if ~strcmp(myFiles(data_index).name,".") && ~strcmp(myFiles(data_index).name,"..") 
            file = readtable(strcat(myDir,'/',myFiles(data_index).name));
                % assume gridding is same for all
            if data_index == 3
                % just for first coil ...
                datasize = size(file,1)^(1/3);
                Ex_allcoils = zeros(30, int8(datasize),int8(datasize),int8(datasize));
                Ey_allcoils = Ex_allcoils;
                Ez_allcoils = Ex_allcoils;  

                X = file.Var1;
                X = X + abs(min(X));
                X = int8((X / max(X))*(datasize-1))+1;
                Y = file.Var2;
                Y = Y + abs(min(Y));
                Y = int8((Y / max(Y))*(datasize-1))+1;
                Z = file.Var3;
                Z = Z + abs(min(Z));
                Z = int8((Z / max(Z))*(datasize-1))+1;
            end
            % data-agnostic, E or H
            Ex = file.Var4 + 1j * file.Var5;
            Ey = file.Var6 + 1j * file.Var7;
            Ez = file.Var8 + 1j * file.Var9;
            
            Ex_reshape = zeros(int8(datasize), int8(datasize),int8(datasize));
            Ey_reshape = Ex_reshape;
            Ez_reshape = Ex_reshape;

            for location_index = 1:size(file,1)
                Ex_reshape(X(location_index,:), Y(location_index,:),Z(location_index,:)) = Ex(location_index,:);
                Ey_reshape(X(location_index,:), Y(location_index,:),Z(location_index,:)) = Ey(location_index,:);
                Ez_reshape(X(location_index,:), Y(location_index,:),Z(location_index,:)) = Ez(location_index,:);
            end
            Ex_allcoils(data_index-2, :,:,:) = Ex_reshape;
            Ey_allcoils(data_index-2, :,:,:) = Ey_reshape;
            Ez_allcoils(data_index-2, :,:,:) = Ez_reshape;
        end
    end
end   

%% do B-field
datadirs = {'head_1_b1field', 'head_2_b1field'};
for dataset_index = 1:length(datadirs)
    % load all of the data in the directory
    myDir = strcat(pwd,'/data/',datadirs{dataset_index}); %gets directory
    myFiles = dir(fullfile(myDir,'*')); %gets all wav files in struct
    
    for data_index = 1:size(myFiles,1)
        data_index
        if ~strcmp(myFiles(data_index).name,".") && ~strcmp(myFiles(data_index).name,"..") 
            file = readtable(strcat(myDir,'/',myFiles(data_index).name));
                % assume gridding is same for all
            if data_index == 3
                % just for first coil ...
                datasize = size(file,1)^(1/3);
                Hx_allcoils = zeros(30, int8(datasize),int8(datasize),int8(datasize));
                Hy_allcoils = Hx_allcoils;
                Hz_allcoils = Hx_allcoils;  

                X = file.Var1;
                X = X + abs(min(X));
                X = int8((X / max(X))*(datasize-1))+1;
                Y = file.Var2;
                Y = Y + abs(min(Y));
                Y = int8((Y / max(Y))*(datasize-1))+1;
                Z = file.Var3;
                Z = Z + abs(min(Z));
                Z = int8((Z / max(Z))*(datasize-1))+1;
            end
            % data-agnostic, E or H
            Hx = file.Var4 + 1j * file.Var5;
            Hy = file.Var6 + 1j * file.Var7;
            Hz = file.Var8 + 1j * file.Var9;
            
            Hx_reshape = zeros(int8(datasize), int8(datasize),int8(datasize));
            Hy_reshape = Hx_reshape;
            Hz_reshape = Hx_reshape;

            for location_index = 1:size(file,1)
                Hx_reshape(X(location_index,:), Y(location_index,:),Z(location_index,:)) = Hx(location_index,:);
                Hy_reshape(X(location_index,:), Y(location_index,:),Z(location_index,:)) = Hy(location_index,:);
                Hz_reshape(X(location_index,:), Y(location_index,:),Z(location_index,:)) = Hz(location_index,:);
            end
            Hx_allcoils(data_index-2, :,:,:) = Hx_reshape;
            Hy_allcoils(data_index-2, :,:,:) = Hy_reshape;
            Hz_allcoils(data_index-2, :,:,:) = Hz_reshape;
        end
    end
    save(strcat(datadirs{dataset_index},'.mat'), 'Hx_allcoils',)
end   

%% Dummy script to generate mask

m = zeros(61,61,61);
for coil_index = 1:15
    sum_coil = squeeze(sqrt(Ex_allcoils(coil_index,:,:,:).^2 + Ey_allcoils(coil_index,:,:,:).^2+Ez_allcoils(coil_index,:,:,:).^2));
    sum_coil(isnan(sum_coil)) = 0;
    m = m + sum_coil;
end
for coil_index = 25:30
    sum_coil = squeeze(sqrt(Ex_allcoils(coil_index,:,:,:).^2 + Ey_allcoils(coil_index,:,:,:).^2+Ez_allcoils(coil_index,:,:,:).^2));
    sum_coil(isnan(sum_coil)) = 0;
    m = m + sum_coil;
end

% segmentation
J = zeros(size(m));
J_nasal_left = J;
J_nasal_right = J;
nasal_thresh = 35;
for slice_index = 1:55
    if slice_index < 6
        T = 35;
    elseif slice_index > 52
        T = 30;
    else
        T = 55;
    end

    J(:,:,slice_index) = regiongrowing(abs(m(:,:,slice_index)),30,25,T);
    if slice_index > 21 && slice_index < 30
        J_nasal_right(:,:,slice_index) = regiongrowing(abs(m(:,:,slice_index)),46,33,nasal_thresh);
        J_nasal_left(:,:,slice_index) = regiongrowing(abs(m(:,:,slice_index)),46,29,nasal_thresh);
    end
    
    % binary operations to fill holes
    J(:,:,slice_index) = imfill(J(:,:,slice_index),"holes");
end
J = J .* ~J_nasal_left .* ~J_nasal_right;



