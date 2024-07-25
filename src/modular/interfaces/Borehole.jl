"""
    abstract type Borehole

Interface for boreholes.

Required functions:
- `get_H(::Borehole)`: Return the length of the borehole.
- `get_D(::Borehole)`: Return the burial depth of the borehole.
- `get_rb(::Borehole)`: Return the radius of the borehole.
"""
abstract type Borehole end

@required Borehole begin
    get_H(::Borehole)
    get_D(::Borehole)
    get_rb(::Borehole)
end
