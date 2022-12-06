function [coil_index_mat] = build_array(number_coils, number_rows, number_cols)
    number_coils = 30; number_rows = 3; number_cols = 10;

    % build the indexing matrix - matches the diagram
    coil_index_mat = [1:1:number_coils];
    coil_index_mat = reshape(coil_index_mat,[number_rows, number_cols]);
    coil_index_mat(3,:) = circshift(coil_index_mat(3,:),1);
    coil_index_mat = flipud(coil_index_mat);

end

