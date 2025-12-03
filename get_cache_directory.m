function cache_directory = get_cache_directory(file_prefix_outputs)
cache_directory = strrep(file_prefix_outputs,'preprocessed_CAP_data','mat_cache_of_raw_data');

if ~exist(cache_directory,'dir')
    mkdir(cache_directory)
end