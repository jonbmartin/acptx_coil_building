function [new_coil_grouping] = acptx_indexing_to_xin_indexing(coil_grouping,lookup_table)
    new_coil_grouping = zeros(size(coil_grouping));
    for ii = 1:size(coil_grouping,1)
        for jj = 1:size(coil_grouping,2)
            if coil_grouping(ii,jj) < 1
                new_coil_grouping(ii,jj) = -1;
            else
                new_coil_grouping(ii,jj) = lookup_table(:,coil_grouping(ii,jj));
            
            end
        end
    end

end

