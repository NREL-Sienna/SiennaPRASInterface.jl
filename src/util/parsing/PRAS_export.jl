function export_pras_system(pras_sys::PRAS.SystemModel,export_location::Union{Nothing, String})
    if (export_location !== nothing)
        if ~(isprasfile(export_location))
            error("PRAS System export location should be a .pras file. $(export_location) is not a valid location.")
        else
            PRAS.savemodel(pras_sys,export_location, string_length =100, verbose = true, compression_level = 9)
            @info "PRAS System exported can be found here : $(export_location)"
        end
    end
end