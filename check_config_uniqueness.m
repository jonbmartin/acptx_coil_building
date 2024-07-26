%groupings of coils. 17 is duplicated in the group with 2 to keep dim
%consistent 
coil_groupings = [1, 2, 3, 30; 4, 5, 6, 7; 8, 9, 10, 11; 12, 13, 14, 15;...
    16, 17, 17, 17; 18, 19, 20, 21; 22, 23, 24, 25; 26,27 28, 29]; 
custom_colormap = coil_colormap(30, coil_groupings);

% build the indexing matrix - matches the diagram
coil_index_mat = build_array(30, 3, 10);

% example: shift 3 and flip vertical

% operations 10 shifts
% vertical flip
% horizontal flip
% both flips 
% AND all rotated 18 dg (separate dataset) 
truth_v = {true, false};
all_coil_combos = {};
combo_index =1
for shift_index = 1:10
    for vert_index = 1:2
        for hor_index = 1:2
            all_coil_combos{combo_index,1} = output_new_config(coil_index_mat, shift_index, truth_v{vert_index}, truth_v{hor_index});
            combo_index = combo_index + 1;
        end
    end
end

% check for uniqueness - expect to return 40 because ii == jj
number_of_duplicates = 0;
for ii = 1:combo_index-1
    for jj = 1:combo_index-1
        number_of_duplicates = number_of_duplicates + isequal(all_coil_combos{ii,:}, all_coil_combos{jj,:});
    end
end
number_of_duplicates = number_of_duplicates - size(all_coil_combos,1); % remove self-comparison


