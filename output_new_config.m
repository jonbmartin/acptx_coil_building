function [coil_index_mat, new_coil_groupings] = output_new_config(coil_index_mat, shift, flip_vertical, flip_horizontal)
%OUTPUT_NEW_CONFIG given a request to move or flip, applies operation to
%coil matrix. Inputs are horizontal shift (positive = right), and boolean
%options to flip the array up or down

    load coil_grouping_ids.mat;

    old_coil_index_mat = coil_index_mat;

    % initial groupings. A 0 is an empty slot (meaning less than 4 coils in a
    % group)
%     coil_groupings = [1, 2, 3, 30; 4, 5, 6, 7; 8, 9, 10, 11; 12, 13, 14, 15;...
%         16, 17, 0, 0; 18, 19, 20, 21; 22, 23, 24, 25; 26,27 28, 29]; 

    % so the coil GROUPINGS stay fixed but the underlying coil index mat is
    % moved around
    if flip_vertical
        coil_index_mat = flipud(coil_index_mat);
    end
    if flip_horizontal
        %coil_index_mat = fliplr(coil_index_mat); % THIS IS NOT CORRECT.
        %MUST SHIFT TO CORRECT
        middle_row = circshift(coil_index_mat(2,:),1);
        coil_index_mat(2,:) = middle_row;
        coil_index_mat = fliplr(coil_index_mat);

    end
    coil_index_mat = circshift(coil_index_mat, shift, 2);
    
    % get the indices in row-sorted groups
    paired_index = [coil_groupings_ids(:), coil_index_mat(:)];
    new_coil_groupings = zeros(8,4);
    uniqueValues = unique(paired_index(:, 1));
    for k = 1 : length(uniqueValues)
        thisValue = uniqueValues(k);
        rowsWithThisValue = paired_index(:, 1) == thisValue;
        subArray = paired_index(rowsWithThisValue, 2);
        % Now sort into 8 x 4
        if length(subArray) == 2
            subArray = [subArray; -1; -1];
        end
        new_coil_groupings(k,:) = subArray;
    end



        


end

