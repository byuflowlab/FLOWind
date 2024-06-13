abstract type AbstractCostModel end

"""
    Levelized(TCC, BOS, FC, FCR, OpEx)

Container for parameters related to the Levelized Cost of Energy model (NREL 2016 Cost of Wind Energy)

# Arguments
# Arguments
- `TCC::Float`: Turbine Capital Cost not including the tower module
- `BOS::Float`: Balance of System (Costs outside of turbine i.e. operation and maintenance)
- `FC::Float`: Financial Costs including construction and contingency
- `FCR::Float`: Fixed Charge Rate
- `OpEx::Float`: Operational Expenditures
"""
struct Levelized{TF} <: AbstractCostModel
    TCC::TF
    BOS::TF
    FC::TF
    FCR::TF
    OpEx::TF
end
Levelized() = Levelized(776.0, 326.0, 120.0, .0655, 43.0) # Default values taken from NREL 2019 Cost of Wind Energy


"""
    cost_of_energy(rotor_diameter, hub_height, AEP, Cost::Levelized)

Calculates the LCOE using the same numbers as NREL's FLORIS Model

# Arguments
- `rotor_diameter::array`: Vector of Rotor Diameters for the Turbines
- `hub_height::array`: Vector of Hub Heights for the Turbines
- `rated_power::array`: Vector of rated powers for the Turbines in kW
- `AEP::Float`: Annual Energy Production in MWh/MW/year
- `OpEx::AbstractCostParameter`: KW of the Farm
"""

function cost_of_energy(rotor_diameter, hub_height, rated_power, AEP, Cost::Levelized)
    
    nturbines = length(rotor_diameter)
    # Taken from 2019 Cost of Wind Energy Review
    TCC = Cost.TCC
    BOS = Cost.BOS
    FC = Cost.FC
    FCR = Cost.FCR
    OpEx = Cost.OpEx

    PlantKW = sum(rated_power) # Combines total Wind Farm Capacity

    Mass = 0.0
    for i = 1:nturbines
        swept_area = pi*(rotor_diameter[i]/2)^2
        Mass += .2694*hub_height[i]*swept_area + 1779.3
    end
    # Adding the mass of the turbines to the Total Capitol Cost
    TCC =  TCC + 3.08*Mass/PlantKW

    # Uses parameters in COE function from eq 1 in 2016 Cost of Wind Energy Review
    LCOE = ((TCC+BOS+FC)*FCR + OpEx)/(AEP/1000)

    return LCOE # LCOE is in units of $/MWh
end
