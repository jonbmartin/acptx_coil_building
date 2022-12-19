%groupings of coils. 17 is duplicated in the group with 2 to keep dim
%consistent 
coil_groupings = [1, 2, 3, 30; 4, 5, 6, 7; 8, 9, 10, 11; 12, 13, 14, 15;...
    16, 17, 17, 17; 18, 19, 20, 21; 22, 23, 24, 25; 26,27 28, 29]; 
custom_colormap = coil_colormap(30, coil_groupings);

% build the indexing matrix - matches the diagram
coil_index_mat = build_array(30, 3, 10);

heatmap(coil_index_mat,FontSize=12);
colormap(custom_colormap)
colorbar off
title('default coil array')

% example: shift 3 and flip vertical
shift = 3; vertical_flip = true; horizontal_flip = false;
coil_index_mat = output_new_config(coil_index_mat, shift, vertical_flip, horizontal_flip);

figure,
heatmap(coil_index_mat,FontSize=12);
colormap(custom_colormap);
colorbar off
title('shifted coil array')
