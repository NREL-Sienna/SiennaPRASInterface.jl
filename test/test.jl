#######################################################
# Surya
# NREL
# January 2021
# SIIP --> PRAS Linkage Module Test
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
pras_system_1 = PowerSystems2PRAS.make_pras_system(sys,system_model="Single-Node",aggregation="Area",period_of_interest=system_period_of_interest,outage_flag=false);
pros_system_1 = PowerSystems2PRAS.make_pras_system(sys_day_ahead,sys_real_time,look_ahead = 3,system_model="Single-Node",aggregation="Area",period_of_interest=system_period_of_interest,outage_flag=true);

### PRAS Util Functions Test  
system_period_of_interest_new = range(10,length = 100);
start_timestamp = pras_system_1.timestamps[10];
stop_timestamp = pras_system_1.timestamps[109];
pras_system_new_1 = PowerSystems2PRAS.reduce_pras_system_size(pras_system_1,system_period_of_interest_new);
pras_system_new_2 = PowerSystems2PRAS.reduce_pras_system_size(pras_system_1,start_timestamp,stop_timestamp);
### PRAS Util Functions Test

pras_system_2 =  PowerSystems2PRAS.make_pras_system(sys,system_model="Zonal",aggregation="Area",period_of_interest=system_period_of_interest,outage_flag=true);
@time pros_system_2 = PowerSystems2PRAS.make_pras_system(sys_day_ahead,sys_real_time,look_ahead = 3,system_model="Zonal",aggregation="Area",period_of_interest=system_period_of_interest,outage_flag=true)
### PRAS Util Functions Test #TODO : Collapse PRAS System when extracted regions are islands
system_period_of_interest_new = range(10,length = 100);
start_timestamp = pras_system_2.timestamps[10];
stop_timestamp = pras_system_2.timestamps[109];
pras_system_2_new_1 = PowerSystems2PRAS.collapse_pras_regions(pras_system_2,period_of_interest=system_period_of_interest_new);
pras_system_2_new_2 = PowerSystems2PRAS.collapse_pras_regions(pras_system_2,start_timestamp=start_timestamp,stop_timestamp=stop_timestamp);
regions_set_1 = [1,3]
regions_set_2 = ["1","3"]
regions_set_3 = [1]
pras_system_2_new_3 = PowerSystems2PRAS.extract_regions(pras_system_2,regions_set_1,period_of_interest=system_period_of_interest_new);
pras_system_2_new_4 = PowerSystems2PRAS.extract_regions(pras_system_2,regions_set_2,period_of_interest=system_period_of_interest_new);
pras_system_2_new_5 = PowerSystems2PRAS.extract_regions(pras_system_2,regions_set_1,start_timestamp=start_timestamp,stop_timestamp=stop_timestamp);
pras_system_2_new_6 = PowerSystems2PRAS.extract_regions(pras_system_2,regions_set_2,start_timestamp=start_timestamp,stop_timestamp=stop_timestamp);

pras_system_2_new_7 = PowerSystems2PRAS.extract_regions(pras_system_2,regions_set_3,system_period_of_interest_new);
### PRAS Util Functions Test

pras_system_3 = PowerSystems2PRAS.make_pras_system(sys,system_model="Single-Node",aggregation="LoadZone",period_of_interest=system_period_of_interest,outage_flag=true);
pras_system_4 = PowerSystems2PRAS.make_pras_system(sys,system_model="Zonal",aggregation="LoadZone",period_of_interest=system_period_of_interest,outage_flag=true);

pras_output_file_location = "/Users/sdhulipa/Desktop/temp/temp-20.pras";
PowerSystems2PRAS.export_pras_system(pras_system_2_new_4,pras_output_file_location);


location = "/Users/sdhulipa/Desktop/OneDrive - NREL/Julia/PowerSystems2PRAS.jl";
PowerSystems2PRAS.generate_outage_profile(pras_system_1,sys, num_runs = 1000,num_scenarios=2,location = location)
PowerSystems2PRAS.generate_csv_outage_profile(pras_system_1, num_runs = 100,num_scenarios=2)

#######################################################
# Testing the module
#######################################################
system_period_of_interest = range(100,length = 500);
pras_system_1 = PowerSystems2PRAS.make_pras_system(sys,"Single-Node","Area",system_period_of_interest);
# @test isa(pras_system_1,PRAS.SystemModel)
pras_system_2 = PowerSystems2PRAS.make_pras_system(sys,"Zonal","Area",system_period_of_interest);
pras_system_3 = PowerSystems2PRAS.make_pras_system(sys,"Single-Node","LoadZone",system_period_of_interest);
pras_system_4 = PowerSystems2PRAS.make_pras_system(sys,"Zonal","LoadZone",system_period_of_interest);

pras_output_file_location = "/Users/sdhulipa/Desktop/temp/temp-4.pras";
PowerSystems2PRAS.export_pras_system(pras_system_1,pras_output_file_location);

# Serialization of HybridSystem is not currently supported
location = "/Users/sdhulipa/Desktop/Julia/PowerSystems2PRAS.jl";
PowerSystems2PRAS.generate_outage_profile(pras_system_1,sys_new, num_runs = 1000,num_scenarios=2,location = location)

pras_system_5 = PowerSystems2PRAS.make_pras_system(sys,"Zonal","Area",100);
PowerSystems2PRAS.generate_outage_profile(pras_system_5,sys, num_runs = 1000,num_scenarios=2,location = location)

#######################################################
# ERCOT System
#######################################################
ercot_data_location = "/Users/sdhulipa/Desktop/temp/Texas RT Bug/texas_data"
sys_DA = PSY.System(joinpath(ercot_data_location,"DA_sys.json"));
sys_RT = PSY.System(joinpath(ercot_data_location,"RT_sys_hourly.json"));
#sys_RT = PSY.System(joinpath(ercot_data_location,"RT_sys.json"));
days_of_interest = range(1,length=364);
look_ahead=35;
@time ercot_pros_system = PowerSystems2PRAS.make_ercot_pras_system(sys_DA,sys_RT,look_ahead=look_ahead ,system_model ="Zonal",aggregation = "Area",
                          days_of_interest =days_of_interest, outage_flag=false,load_DA_look_ahead_weight=1.0,solar_DA_look_ahead_weight = 1.0);

ercot_rt_pros_system = PowerSystems2PRAS.make_RT_ercot_pras_system(sys_DA,sys_RT,look_ahead = look_ahead,system_model= "Zonal",aggregation= "Area",days_of_interest= days_of_interest,
outage_flag=false,lump_wind_gens=true,load_DA_look_ahead_weight = 1.0,solar_DA_look_ahead_weight= 1.0);

# Workflow to add generator availability data to DA and RT PSY Systems (Two Stage SIIP Simulations)
ercot_data_location = "/Users/sdhulipa/Desktop/temp/Texas RT Bug/texas_data"
sys_DA = PSY.System(joinpath(ercot_data_location,"DA_sys.json"));
sys_RT = PSY.System(joinpath(ercot_data_location,"RT_sys.json"));
days_of_interest = range(1,length=3);
look_ahead=35;
@time ercot_rt_pros_system = PowerSystems2PRAS.make_RT_ercot_pras_system(sys_DA,sys_RT,look_ahead = look_ahead,system_model= "Zonal",aggregation= "Area",days_of_interest= days_of_interest,
outage_flag=false,lump_wind_gens=true,load_DA_look_ahead_weight = 1.0,solar_DA_look_ahead_weight= 1.0);

#######################################################
# Saving and Loading the PrOS System
#######################################################
ercot_pras_file_location = "/Users/sdhulipa/Desktop/OneDrive - NREL/Julia/PowerSystems2PRAS.jl/data/ERCOT-RT-PrOS-365-35.jld"

PowerSystems2PRAS.export_pros_system(ercot_rt_pros_system,ercot_pras_file_location,pros_sys_name="ercot_rt_365_35")
import PrOS
ercot_rt_pros_system  = PowerSystems2PRAS.load_pros_system(ercot_pras_file_location,pros_sys_name="ercot_rt_365_35");


num_runs = 10
result = PrOS.assess(ercot_rt_pros_system, PrOS.EconomicDispatch(samples = num_runs),PrOS.UnservedSummary(), PrOS.UnservedSamples(), PrOS.UnservedFull(), PrOS.Availability(), PrOS.StorageDispatch(), PrOS.Transmission()) 
outage_csv_location = PowerSystems2PRAS.extract_csv_outage_profile(result,ercot_rt_pros_system,num_scenarios = 4)
sys_RT = PSY.System(joinpath(ercot_data_location,"RT_sys_hourly.json"));
new_sys_DA, new_sys_RT = PowerSystems2PRAS.add_csv_time_series!(sys_DA,sys_RT,outage_csv_location,add_scenario = 1);
# Workflow to add generator availability data to DA System (Single Stage SIIP Simulations)
new_sys_DA = PowerSystems2PRAS.add_csv_time_series_single_stage!(sys_DA,outage_csv_location,add_scenario = 3);
#######################################################
# Saving and Loading the PrOS System
#######################################################
ercot_pras_file_location = "/Users/sdhulipa/Desktop/OneDrive - NREL/Julia/PowerSystems2PRAS.jl/data/ERCOT-RAW-PRAS-ED-365-3.jld"

PowerSystems2PRAS.export_pros_system(ercot_pros_system,ercot_pras_file_location,pros_sys_name="ercot_raw_365_3")
import PrOS
ercot_365_3_pros_system = PowerSystems2PRAS.load_pros_system(ercot_pras_file_location,pros_sys_name="ercot_365_3");

# cc_restrictions and outage CSV's for testing
Location = "/Users/sdhulipa/Desktop/temp/Texas RT Bug/texas_data"
Location_1 = joinpath(Location,"cc_restrictions.json")
cc_restrictions = JSON.parsefile(Location_1);

Location_2 = "/Users/sdhulipa/Desktop/OneDrive - NREL/Julia/PowerSystems2PRAS.jl/src/main/descriptors"; 
Location_3 = joinpath(Location_2,"outage-rates-ERCOT-modified.csv")
df_outage = DataFrames.DataFrame(CSV.File(Location_3));


