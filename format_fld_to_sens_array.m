function [] = format_fld_to_sens_array(stringToFile)
%FORMAT_FLD_TO_ Summary of this function goes here
%   Detailed explanation goes here
    load(stringToFile);
    X = chExEyEz1{:,1};
    Y = chExEyEz1{:,2};
    Z = chExEyEz1{:,3};
    Ex = chExEyEz1{:,4} + 1j * chExEyEz1{:,5};
    Ey = chExEyEz1{:,4} + 1j * chExEyEz1{:,5};
    Ez = chExEyEz1{:,4} + 1j * chExEyEz1{:,5};


end

