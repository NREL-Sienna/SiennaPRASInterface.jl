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
rts_dir = download("https://github.com/scdhulipala/RTS-GMLC", "master", joinpath(pwd(),"data"))
rts_src_dir = joinpath(rts_dir, "RTS_Data", "SourceData");
rts_siip_dir = joinpath(rts_dir, "RTS_Data", "FormattedData", "SIIP");

rawsys = PSY.PowerSystemTableData(
    rts_src_dir,
    100.0,
    joinpath(rts_siip_dir, "user_descriptors.yaml"),
    timeseries_metadata_file = joinpath(rts_siip_dir, "timeseries_pointers.json"),
    generator_mapping_file = joinpath(rts_siip_dir, "generator_mapping.yaml"),
);

sys = PSY.System(rawsys; time_series_resolution = Dates.Hour(1));
#######################################################
# PrOS Stuff
#######################################################
sys_day_ahead = PSY.System(rawsys; time_series_resolution = Dates.Hour(1));
sys_real_time = PSY.System(rawsys; time_series_resolution = Dates.Minute(5));
#######################################################
# Testing the module
# You have to load the module before
#######################################################
system_period_of_interest = range(1,length = 8782);
pras_system_1 = SIIP2PRAS.make_pras_system(sys,system_model="Single-Node",aggregation="Area",period_of_interest=system_period_of_interest,outage_flag=false);
pras_system_2 =  SIIP2PRAS.make_pras_system(sys,system_model="Zonal",aggregation="Area",period_of_interest=system_period_of_interest,outage_flag=true);

pras_system_3 = SIIP2PRAS.make_pras_system(sys,system_model="Single-Node",aggregation="LoadZone",period_of_interest=system_period_of_interest,outage_flag=true);
pras_system_4 = SIIP2PRAS.make_pras_system(sys,system_model="Zonal",aggregation="LoadZone",period_of_interest=system_period_of_interest,outage_flag=true);


location = "/Users/sdhulipa/Desktop/OneDrive - NREL/Julia/PowerSystems2PRAS.jl";
SIIP2PRAS.generate_outage_profile(pras_system_1,sys, num_runs = 1000,num_scenarios=2,location = location)
SIIP2PRAS.generate_csv_outage_profile(pras_system_1, num_runs = 100,num_scenarios=2)

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



