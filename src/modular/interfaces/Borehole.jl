"""
    abstract type Borehole

Common interface for boreholes
"""
abstract type Borehole end

#  Get the length of the borehole
function get_H(::Borehole) end
# Get the burial depth of the borehole
function get_D(::Borehole) end
# Get the segment length of the borehole
function get_h(::Borehole) end
# Get the radius of the borehole
function get_rb(::Borehole) end