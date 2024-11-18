# TODO: OUTAGE_INFO FILE: this should probably be an artifact 
"""
    DEFAULT outage data which is used when outage_flag is set to FALSE

Based on ERCOT historical data
"""
const OUTAGE_INFO_FILE =
    joinpath(@__DIR__, "descriptors", "outage-rates-ERCOT-modified.csv")

"""
    Supported DC branch Types
"""
const HVDCLineTypes = Union{PSY.TwoTerminalHVDCLine, PSY.TwoTerminalVSCDCLine}

"""
    Supported Transformer Types
"""
const TransformerTypes =
    [PSY.TapTransformer, PSY.Transformer2W, PSY.PhaseShiftingTransformer]
