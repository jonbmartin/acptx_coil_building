function [colormap] = coil_colormap(number_coils, grouping)
%COIL_COLORMAP Summary of this function goes here
%   Detailed explanation goes here
    colormap = zeros(number_coils, 3);
    number_of_groups = size(grouping,1);
    for index = 1:number_of_groups
        colormap(grouping(index, :),:) = index/number_of_groups;
        colormap(grouping(index, :),:) = index/number_of_groups;
    end
end

