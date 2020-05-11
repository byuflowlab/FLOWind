abstract type AbstractPowerModel end

struct PowerModelConstantCp{TF} <: AbstractPowerModel
    cp::TF
end

struct PowerModelCpPoints{ATF} <: AbstractPowerModel
    vel_points::ATF
    cp_points::ATF
end

struct PowerModelPowerPoints{ATF} <: AbstractPowerModel
    vel_points::ATF
    power_points::ATF
end

struct PowerModelPowerCurveCubic{} <: AbstractPowerModel
end

"""
    calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity)

Calculate the power for a wind turbine based on standard theory for region 2

# Arguments
- `generator_efficiency::Float`: Efficiency of the turbine generator
- `air_density::Float`: Air density
- `rotor_area::Float`: Rotor-swept area of the wind turbine
- `cp::Float`: Power coefficient of the wind turbine
- `wt_velocity::Float`: Inflow velocity to the wind turbine
"""
function calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity)
    power = generator_efficiency*(0.5*air_density*rotor_area*cp*wt_velocity^3)
    return power
end

"""
    calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, power_model)

Calculate the power for a wind turbine based on standard theory for region 2 using a constant cp

# Arguments
- `generator_efficiency::Float`: Efficiency of the turbine generator
- `air_density::Float`: Air density
- `rotor_area::Float`: Rotor-swept area of the wind turbine
- `wt_velocity::Float`: Inflow velocity to the wind turbine
- `cut_in_speed::Float`: cut in speed of the wind turbine
- `rated_speed::Float`: rated speed of the wind turbine
- `cut_out_speed::Float`: cut out speed of the wind turbine
- `rated_power::Float`: rated power of the wind turbine
- `power_model::PowerModelConstantCp`: Struct containing the cp value to be used in region 2
"""
function calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model::PowerModelConstantCp)

    if wt_velocity < cut_in_speed
        power = 0.0
    elseif wt_velocity < rated_speed
        # extract cp_value
        cp = power_model.cp
        power = calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity)
        if power > rated_power
            power = rated_power
        end
    elseif wt_velocity < cut_out_speed
        power = rated_power
    elseif wt_velocity > cut_out_speed
        power = 0.0
    end

    return power

end

"""
    calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model)

Calculate the power for a wind turbine based on a cp curve with linear interpolation

# Arguments
- `generator_efficiency::Float`: Efficiency of the turbine generator
- `air_density::Float`: Air density
- `rotor_area::Float`: Rotor-swept area of the wind turbine
- `wt_velocity::Float`: Inflow velocity to the wind turbine
- `cut_in_speed::Float`: cut in speed of the wind turbine
- `rated_speed::Float`: rated speed of the wind turbine
- `cut_out_speed::Float`: cut out speed of the wind turbine
- `rated_power::Float`: rated power of the wind turbine
- `power_model::PowerModelCpPoints`: Struct containing the velocity and cp values defining the cp curve
"""
function calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model::PowerModelCpPoints)

    # use specs if inflow wind speed is less than the wind speeds provided in the power curve
    if wt_velocity < power_model.vel_points[1]
        
        # calculated wind turbine power
        if wt_velocity < cut_in_speed
            power = 0.0
        elseif wt_velocity < rated_speed
            # use cp value corresponding to lowest provided velocity point
            cp = power_model.cp_points[1]
            # calculate power
            power = calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity) 
        end

    # use cp points where provided
    elseif wt_velocity < power_model.vel_points[end]

        # estimate cp_value using linear interpolation
        cp = linear(power_model.vel_points, power_model.cp_points, wt_velocity)

        # calculate power
        power = calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity) 
        
    # use specs if above vel_points max
    else

        if wt_velocity <= cut_out_speed
            # use cp value corresponding to highest provided velocity point
            cp = power_model.cp_points[end]
            # calculate power
            power = calculate_power_from_cp(generator_efficiency, air_density, rotor_area, cp, wt_velocity) 
        elseif wt_velocity > cut_out_speed
            power = 0.0
        end

    end

    return power

end

"""
    calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model)

Calculate the power for a wind turbine based on a pre-determined power curve with linear
    interpolation

# Arguments
- `generator_efficiency::Float`: Efficiency of the turbine generator
- `air_density::Float`: Air density
- `rotor_area::Float`: Rotor-swept area of the wind turbine
- `wt_velocity::Float`: Inflow velocity to the wind turbine
- `cut_in_speed::Float`: cut in speed of the wind turbine
- `rated_speed::Float`: rated speed of the wind turbine
- `cut_out_speed::Float`: cut out speed of the wind turbine
- `rated_power::Float`: rated power of the wind turbine
- `power_model::PowerModelPowerPoints`: Struct containing the velocity and power values
    defining the power curve
"""
function calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model::PowerModelPowerPoints)

    # use specs if inflow wind speed is less than the wind speeds provided in the power curve
    if wt_velocity < power_model.vel_points[1]
        
        # calculated wind turbine power
        if wt_velocity < cut_in_speed
            power = 0.0
        elseif wt_velocity < rated_speed
            # use power value corresponding to lowest provided velocity point
            power = linear([cut_in_speed, power_model.vel_points[1]], [0.0, power_model.power_points[1]], wt_velocity)
        end

    # use power points where provided
    elseif wt_velocity < power_model.vel_points[end]

        # calculate power
        power = linear(power_model.vel_points, power_model.power_points, wt_velocity)
        
    # use specs if above vel_points max
    else

        if wt_velocity <= cut_out_speed
            # use power corresponding to highest wind speed provided
            power = power_model.power_points[end]
        elseif wt_velocity > cut_out_speed
            power = 0.0
        end

    end

    return power

end

"""
    calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model)

Calculates wind turbine power using a cubic estimation based on turbine specifications
    as defined in https://github.com/byuflowlab/iea37-wflo-casestudies/blob/master/cs3-4/iea37-cs3-announcement.pdf

# Arguments
- `generator_efficiency::Float`: Efficiency of the turbine generator
- `air_density::Float`: Air density
- `rotor_area::Float`: Rotor-swept area of the wind turbine
- `wt_velocity::Float`: Inflow velocity to the wind turbine
- `cut_in_speed::Float`: cut in speed of the wind turbine
- `rated_speed::Float`: rated speed of the wind turbine
- `cut_out_speed::Float`: cut out speed of the wind turbine
- `rated_power::Float`: rated power of the wind turbine
- `power_model::PowerModelPowerCurveCubic`: Empty struct
"""
function calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model::PowerModelPowerCurveCubic)
        
    
    if wt_velocity < cut_in_speed
        power = 0.0
    elseif wt_velocity < rated_speed
        power = rated_power*((wt_velocity - cut_in_speed)/(rated_speed - cut_in_speed))^3
    elseif wt_velocity < cut_out_speed
        power = rated_power
    elseif wt_velocity > cut_out_speed
        power = 0.0
    end

    return power

end

"""
    calculate_turbine_power(turbine_id, turbine_definition::TurbineDefinition,
        farmstate::SingleWindFarmState, wind_model::AbstractWindResourceModel)

Calculate the power for a wind turbine based on a pre-determined power curve with linear
    interpolation

# Arguments
- `turbine_id::Int`: Efficiency of the turbine generator
- `turbine_definition::TurbineDefinition`: Struct containing the relevant wind turbine deffinition
- `farmstate::SingleWindFarmState`: Struct contatining the current wind farm state, including correct inflow velocities
- `wind_model::AbstractWindResourceModel`: Struct defining the wind resource
"""
function calculate_turbine_power(turbine_id, turbine_definition::TurbineDefinition, farmstate::SingleWindFarmState, wind_model::AbstractWindResourceModel)

    # extract turbine design information
    generator_efficiency = turbine_definition.generator_efficiency[1]
    cut_in_speed = turbine_definition.cut_in_speed[1]
    cut_out_speed = turbine_definition.cut_out_speed[1]
    rated_speed = turbine_definition.rated_speed[1]
    rated_power = turbine_definition.rated_power[1]
    rotor_diameter = turbine_definition.rotor_diameter[1]

    # calculated wind turbine rotor-swept area
    rotor_area = pi*(rotor_diameter^2)/4.0

    # extract wind resource information
    air_density = wind_model.air_density

    # extract wind turbine inflow velocity
    wt_velocity = farmstate.turbine_inflow_velcities[turbine_id]

    # extract wind turbine power model
    power_model = turbine_definition.power_model

    wt_power = calculate_power(generator_efficiency, air_density, rotor_area, wt_velocity, cut_in_speed, rated_speed, cut_out_speed, rated_power, power_model)

    return wt_power
end


function turbine_powers_one_direction!(rotor_sample_points_y, rotor_sample_points_z,
    problem_description::AbstractWindFarmProblem; wind_farm_state_id=1)

    windfarm = problem_description.wind_farm
    wind_model = problem_description.wind_resource
    farmstate = problem_description.wind_farm_states[wind_farm_state_id]

    # get number of turbines and rotor sample point
    n_turbines = length(farmstate.turbine_x)

    for d=1:n_turbines
        turbine = windfarm.turbine_definitions[windfarm.turbine_definition_ids[d]]
        wt_power = calculate_turbine_power(d, turbine, farmstate, wind_model)
        farmstate.turbine_generators_powers[d] = wt_power
    end
end


"""
    calculate_aep(model_set::AbstractModelSet, problem_description::AbstractWindFarmProblem;
            rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0])

Calculate wind farm AEP


# Arguments
- `turbine_id::Int`: Efficiency of the turbine generator
- `turbine_definition::TurbineDefinition`: Struct containing the relevant wind turbine deffinition
- `farmstate::SingleWindFarmState`: Struct contatining the current wind farm state, including correct inflow velocities
- `wind_model::AbstractWindResourceModel`: Struct defining the wind resource
"""
function calculate_aep(model_set::AbstractModelSet, problem_description::AbstractWindFarmProblem;
            rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0])

    """NEED TESTS"""
    wind_farm = problem_description.wind_farm
    wind_resource = problem_description.wind_resource
    wind_probabilities = wind_resource.wind_probabilities

    nstates = length(problem_description.wind_farm_states)
    hours_per_year = 365.25*24.0

    state_energy = zeros(nstates)
    for i = 1:nstates
        problem_description.wind_farm_states[i].turbine_x[:],problem_description.wind_farm_states[i].turbine_y[:] =
                rotate_to_wind_direction(wind_farm.turbine_x, wind_farm.turbine_y, wind_resource.wind_directions[i])

        problem_description.wind_farm_states[i].sorted_turbine_index[:] = sortperm(problem_description.wind_farm_states[i].turbine_x)

        turbine_velocities_one_direction!(rotor_sample_points_y, rotor_sample_points_z,
            model_set, problem_description, wind_farm_state_id=i)

        turbine_powers_one_direction!(rotor_sample_points_y, rotor_sample_points_z,
            problem_description, wind_farm_state_id=i)

        state_power = sum(problem_description.wind_farm_states[i].turbine_generators_powers)
        state_energy[i] = state_power*hours_per_year*wind_probabilities[i]
    end

    return sum(state_energy)
end
