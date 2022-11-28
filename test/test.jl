#######################################################
# Surya
# NREL
# August 2022
# SIIP2PRAS Module Test
#######################################################
# Loading the required packages
#######################################################
import PowerSystems
const PSY = PowerSystems
import Dates

#import module here
#######################################################
# System from tabular Data
# ### Fetch Data
# ### The tabular data format relies on a folder containing `*.csv` files and a `user_descriptors.yaml` file
# First, we'll read the tabular data
#######################################################
rts_dir = "/Users/sdhulipa/Desktop/Misc./temp/RTS-GMLC"
rts_src_dir = joinpath(rts_dir, "RTS_Data", "SourceData");
rts_siip_dir = joinpath(rts_dir, "RTS_Data", "FormattedData", "SIIP");
user_descriptrs_yaml_location = joinpath("/Users/sdhulipa/Desktop/OneDrive - NREL/NREL-Github/SIIP2PRAS/test/descriptors/user_descriptors.yaml") # you have to use this yaml if you want to parse FOR and MTTR

rawsys = PSY.PowerSystemTableData(
    rts_src_dir,
    100.0,
    user_descriptrs_yaml_location,
    timeseries_metadata_file = joinpath(rts_siip_dir, "timeseries_pointers.json"),
    generator_mapping_file = joinpath(rts_siip_dir, "generator_mapping.yaml"),
);

sys = PSY.System(rawsys; time_series_resolution = Dates.Hour(1)); # This system should have "outage_probability" and "recovery_probability" in ext of components with outage data.
#######################################################
# Testing the module
# You have to load the module before
#######################################################
system_period_of_interest = range(1,length = 300);
pras_system_1 = SIIP2PRAS.make_pras_system(sys,system_model="Single-Node",aggregation="Area",period_of_interest=system_period_of_interest,outage_flag=false,lump_pv_wind_gens=false,availability_flag=false);
pras_system_2 = SIIP2PRAS.make_pras_system(sys,system_model="Single-Node",aggregation="Area",period_of_interest=system_period_of_interest,outage_flag=false,lump_pv_wind_gens=true);
pras_system_3 =  SIIP2PRAS.make_pras_system(sys,system_model="Zonal",aggregation="Area",period_of_interest=system_period_of_interest,outage_flag=true,lump_pv_wind_gens=false);


location = dirname(pwd());
SIIP2PRAS.generate_outage_profile(pras_system_1,sys, num_runs = 100,num_scenarios=2,location = location)
SIIP2PRAS.generate_csv_outage_profile(pras_system_2, num_runs = 100,num_scenarios=1)

#######################################################
# Testing the module
#######################################################
system_period_of_interest = range(100,length = 500);
pras_system_1 = SIIP2PRAS.make_pras_system(sys,"Single-Node","Area",system_period_of_interest);
# @test isa(pras_system_1,PRAS.SystemModel)
pras_system_2 = SIIP2PRAS.make_pras_system(sys,"Zonal","Area",system_period_of_interest);
pras_system_3 = SIIP2PRAS.make_pras_system(sys,"Single-Node","LoadZone",system_period_of_interest);
pras_system_4 = SIIP2PRAS.make_pras_system(sys,"Zonal","LoadZone",system_period_of_interest);

pras_output_file_location = "/Users/sdhulipa/Desktop/temp/temp-4.pras";
SIIP2PRAS.export_pras_system(pras_system_1,pras_output_file_location);

# Serialization of HybridSystem is not currently supported
location = "/Users/sdhulipa/Desktop/Julia/PowerSystems2PRAS.jl";
SIIP2PRAS.generate_outage_profile(pras_system_1,sys_new, num_runs = 1000,num_scenarios=2,location = location)

pras_system_5 = SIIP2PRAS.make_pras_system(sys,"Zonal","Area",100);
SIIP2PRAS.generate_outage_profile(pras_system_5,sys, num_runs = 1000,num_scenarios=2,location = location)

# Yinong Test System
DA_location = "/Users/sdhulipa/Desktop/OneDrive - NREL/NREL-Github/SIIP2PRAS/data/Test_Data/DA_sys_EMIS_v0811.json" 
sys_DA = PSY.System(DA_location);

outage_csv_location = "/Users/sdhulipa/Desktop/OneDrive - NREL/NREL-Github/SIIP2PRAS/data/Test_Data/ThermalFOR_2012.csv"

pras_system = SIIP2PRAS.make_pras_system(sys_DA,system_model="Single-Node",aggregation="Area",outage_flag=false,lump_pv_wind_gens=false,availability_flag=true,outage_csv_location = outage_csv_location);

pras_system = SIIP2PRAS.make_pras_system(sys_DA,system_model="Single-Node",aggregation="Area",outage_flag=false,lump_pv_wind_gens=false,availability_flag=true);

pras_system = SIIP2PRAS.make_pras_system(sys_DA,system_model="Single-Node",aggregation="Area",outage_flag=true,lump_pv_wind_gens=false,availability_flag=true);