# evaluate nusselt number for flow in pipe for given value of Reynolds and Prandtl number
# currently for Reynolds greater than 10000 we use Dittus-Boelter correlation which require
# knowledge on whether the fluid is heated or cooled. The choosen default is that the fluid is heated.
function evaluate_nusselt(Re,Pr; mode = "heating")
    if Re<2300
        Nu = 3.66
    elseif (Re >= 2300) && (Re<10000)        
        f =  (1.58*log(Re) - 3.28)^(-2)
        Nu = (0.5*f*(Re - 1000.) * Pr)/(1 + 12.7*sqrt(0.5*f) *(Pr^(2/3) - 1))
    else 
        n = mode == "heating" ? 0.3 : 0.4  
        Nu = 0.023 * Re^0.8 * Pr^n
    end
    return Nu
end

# evaluate fluid propert at a given absolute temperature. 
# default mixture is water&Ethanol with 20% concentration of Ethanol
function thermophysical_properties(Tref, fluidname::String)
    μ =  PropsSI("viscosity","T",Tref,"P",101325,fluidname)
    ρ =  PropsSI("D","T", Tref,"P",101325,fluidname)
    cp = PropsSI("C","T",Tref,"P",101325,fluidname)
    k  = PropsSI("conductivity","T",Tref,"P",101325,fluidname)	
    return μ, ρ, cp, k
end

# Evaluate the heat transfer coefficient 
function heat_transfer_coefficient(mf, Tref, borehole::Borehole, fluid)
    T0 = 273.15
    rp = get_rp(borehole)
    μ, ρ, cp, k = thermophysical_properties(fluid, Tref + T0)
    w = mf/(ρ * π * rp^2)
    Re = 2 * ρ * w * rp/ μ
    Pr = μ * cp/(2*rp)
    Nu = evaluate_nusselt(Re, Pr)
    return Nu * k /(2*rp)
end
