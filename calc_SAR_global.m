function [Sglobal] = calc_SAR_global(Cond3dxyz,MassDensity3dxyz, Ex3dxyz, Ey3dxyz, Ez3dxyz, Mask3dxyz)

    %get global SAR matrix
    SAR_factor = Cond3dxyz./MassDensity3dxyz;
    SAR_factor(isnan(SAR_factor)) = 0;
    SAR_factor(isinf(SAR_factor)) = 0;
    Sx = sqrt(0.5)*Ex3dxyz.*repmat(sqrt(SAR_factor),[1 1 1 Nc]);
    Sy = sqrt(0.5)*Ey3dxyz.*repmat(sqrt(SAR_factor),[1 1 1 Nc]);
    Sz = sqrt(0.5)*Ez3dxyz.*repmat(sqrt(SAR_factor),[1 1 1 Nc]);
    Sx = permute(Sx,[4 1 2 3]);Sx = Sx(:,:).';Sx = Sx(logical(Mask3dxyz),:);
    Sy = permute(Sy,[4 1 2 3]);Sy = Sy(:,:).';Sy = Sy(logical(Mask3dxyz),:);
    Sz = permute(Sz,[4 1 2 3]);Sz = Sz(:,:).';Sz = Sz(logical(Mask3dxyz),:);
    Ns = size(Sx,1);
    S = zeros(Nc,Nc,Ns);
    for ii = 1:Ns
        S(:,:,ii) = Sx(ii,:)'*Sx(ii,:) + Sy(ii,:)'*Sy(ii,:) + Sz(ii,:)'*Sz(ii,:);
    end
    Sglobal = mean(S,3);
end

