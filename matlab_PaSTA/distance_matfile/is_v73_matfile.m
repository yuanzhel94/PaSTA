function flag = is_v73_matfile(file)
    flag = false;
    if ~isfile(file)
        return
    end
    try
        info = h5info(file);
        flag = ~isempty(info);  % v7.3 MAT-files use HDF5 format
    catch
        flag = false;
    end
end