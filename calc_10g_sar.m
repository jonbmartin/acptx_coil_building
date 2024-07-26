function [S] = calc_10g_sar(maps)
    % calculate the 10g SAR maps
    Nc = size(maps.Ex_out,4);
    Ngram = 10;
    % note normalize by singular values times Nc divide by sum.. check 

    % want to develop safeguards... 
    
    
    Ex = maps.Ex_out;
    Ey = maps.Ey_out;
    Ez = maps.Ez_out;
    cond3D = maps.cond3D_out;
    
    dens3D = maps.dens3D_out;
    
    SAR_factor = cond3D ./ dens3D;
    SAR_factor(isnan(SAR_factor)) = 0;
    SAR_factor(isinf(SAR_factor)) = 0;
    
    clear *_out;
    Sx = sqrt(0.5)*Ex.*repmat(sqrt(SAR_factor),[1 1 1 Nc]);
    Sy = sqrt(0.5)*Ey.*repmat(sqrt(SAR_factor),[1 1 1 Nc]);
    Sz = sqrt(0.5)*Ez.*repmat(sqrt(SAR_factor),[1 1 1 Nc]);
    mask = sum(abs(Ex) + abs(Ey) + abs(Ez),4).*abs(SAR_factor) > 0;
    Sx = permute(Sx,[4 1 2 3]);Sx = Sx(:,:).';Sx = Sx(mask,:);
    Sy = permute(Sy,[4 1 2 3]);Sy = Sy(:,:).';Sy = Sy(mask,:);
    Sz = permute(Sz,[4 1 2 3]);Sz = Sz(:,:).';Sz = Sz(mask,:);
    
    Ns = size(Sx,1);
    S = zeros(Nc,Nc,Ns);
    for ii = 1:Ns
        S(:,:,ii) = Sx(ii,:)'*Sx(ii,:) + Sy(ii,:)'*Sy(ii,:) + Sz(ii,:)'*Sz(ii,:);
    end
end

