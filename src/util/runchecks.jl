isjson = endswith(".json")
iscsv = endswith(".csv")
isprasfile = endswith(".pras")

"""
    check_file(loc::String)

Check if the file exists and is openable.
"""
function check_file(loc::String)
    io = try
        open(loc)
    catch
        nothing
    end

    if (isnothing(io))
        return nothing, false
    else
        return io, isopen(io)
    end
end

"""
    runchecks(sys_location::String)

Check if the System JSON file is serialized as well as other files required.
"""
function runchecks(sys_location::String)
    if ~(isjson(sys_location))
        error("The System location passed is not a serialized System JSON file")
    end

    sys_dir = dirname(sys_location)
    sys_json_name = basename(sys_location)
    other_sienna_files = [
        "$(first(split(sys_json_name,".json")))_time_series_storage.h5",
        "$(first(split(sys_json_name,".json")))_metadata.json",
    ]
    (root, dirs, files) = first(walkdir(sys_dir))
    if (~all(other_sienna_files .∈ Ref(files)))
        error(
            "All files avaialble to de-serialize the System aren't available in the System JSON folder",
        )
    end

    io, bool = check_file(
        joinpath(sys_dir, "$(first(split(sys_json_name,".json")))_time_series_storage.h5"),
    )
    if (bool)
        close(io)
    else
        error(
            "You don't have access to $(first(split(sys_json_name,".json")))_time_series_storage.h5. Change permissions to this file.",
        )
    end
end
