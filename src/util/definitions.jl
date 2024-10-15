# DEFAULT outage data which is used when outage_flag is set to FALSE
# Baed on ERCOT historical data
const OUTAGE_INFO_FILE =
    joinpath(@__DIR__, "descriptors", "outage-rates-ERCOT-modified.csv") 

# ACBranch Types
const HVDCLineTypes = Union{PSY.TwoTerminalHVDCLine, PSY.TwoTerminalVSCDCLine}
const TransformerTypes = [PSY.TapTransformer, PSY.Transformer2W,PSY.PhaseShiftingTransformer]