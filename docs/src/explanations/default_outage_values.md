# Default values used for outage statistics

When data is not available in a user-supplied `GeometricDistributionForcedOutage`
supplemental attribute, then the outage rates default to a set of defaults defined
in the [Default Outage Rates CSV](https://github.com/NREL-Sienna/SiennaPRASInterface.jl/blob/main/src/util/descriptors/outage-rates-ERCOT-modified.csv) based off of rates in ERCOT.

For any remaining components not captured by the supplemental attributes or CSV,
such as lines, renewables, and storage, the outage rates are 0 0 and will never
fail.