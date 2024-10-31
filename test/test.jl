#######################################################
# Surya
# NREL
# October 2024
# Sienna2PRAS Module Test
#######################################################
# Use PSCB to build RTS-GMLC System and add outage data
# as a Supplmental Attribute
###########################################
import PowerSystems
import PowerSystemCaseBuilder
import CSV
import DataFrames

const PSY = PowerSystems
const PSCB = PowerSystemCaseBuilder

rts_da_sys = PSCB.build_system(PSCB.PSISystems, "RTS_GMLC_DA_sys");
PSY.set_units_base_system!(rts_da_sys, "natural_units")
###########################################
# Parse the gen.csv and add OutageData 
# SupplementalAttribute to components for 
# which we have this data 
###########################################
gen_for_data = CSV.read("descriptors/gen.csv", DataFrames.DataFrame);

for row in DataFrames.eachrow(gen_for_data)
    λ, μ = rate_to_probability(row.FOR, row["MTTR Hr"])
    transition_data = PSY.GeometricDistributionForcedOutage(;
        mean_time_to_recovery=row["MTTR Hr"],
        outage_transition_probability=λ,
    )
    comp = PSY.get_component(PSY.Generator, rts_da_sys, row["GEN UID"])

    if ~(isnothing(comp))
        PSY.add_supplemental_attribute!(rts_da_sys, comp, transition_data)
        @info "Added outage data supplemental attribute to $(row["GEN UID"]) generator"
    else
        @warn "$(row["GEN UID"]) generator doesn't exist in the System."
    end
end

# Make a PRAS System from PSY-4.X System
rts_pras_sys = Sienna2PRAS.make_pras_system(rts_da_sys, PSY.Area);
rts_pras_sys =
    Sienna2PRAS.make_pras_system(rts_da_sys, PSY.Area, lump_region_renewable_gens=true);
rts_pras_sys = Sienna2PRAS.make_pras_system(
    rts_da_sys,
    PSY.Area,
    lump_region_renewable_gens=true,
    availability=false,
);
rts_pras_sys = Sienna2PRAS.make_pras_system(
    rts_da_sys,
    PSY.Area,
    lump_region_renewable_gens=true,
    availability=false,
    export_location="rts.pras",
);
rts_sys_location = "/Users/sdhulipa/Old Mac Backup/Desktop/OneDrive-Backup/NREL-Github/temp/PSCB_test/RTS_4.X.json"
rts_pras_sys = Sienna2PRAS.make_pras_system(
    rts_sys_location,
    PSY.Area,
);
