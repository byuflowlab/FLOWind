abstract type AbstractModelSet end
# using CSV
# using DataFrames

"""
    WindFarmModelSet(wakedeficitmodel, wake_deflection_model, wake_combination_model, local_ti_model)

Container for objects defining models to use in wind farm calculations

# Arguments
- `wake_defiict_model::AbstractWakeDeficitModel`: contains a struct defining the desired wake deficit model
- `wake_deflection_model::AbstractWakeDeflectionModel`: contains a struct defining the desired wake deflection model
- `wake_combination_model::AbstractWakeCombinationModel`: contains a struct defining the desired wake combination model
- `local_ti_model::AbstractTurbulenceIntensityModel`: contains a struct defining the desired turbulence intensity model
"""
struct WindFarmModelSet{DTM,DNM,CM,TIM} <: AbstractModelSet

    wake_deficit_model::DTM
    wake_deflection_model::DNM
    wake_combination_model::CM
    local_ti_model::TIM

end

"""
    point_velocity(loc, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
    wind_resource, model_set::AbstractModelSet;
    wind_farm_state_id=1, downwind_turbine_id=0)

Calculates the wind speed at a given point for a given state

# Arguments
- `loc::Array{TF,3}`: Location of interest
- `turbine_x::Array{TF,nTurbines}`: turbine east-west locations in the state 
    reference frame
- `turbine_y::Array{TF,nTurbines}`: turbine north-south locations in the state 
    reference frame
- `turbine_z::Array{TF,nTurbines}`: turbine base height in the state reference frame
- `turbine_yaw::Array{TF,nTurbines}`: turbine yaw for the given wind direction in 
    radians
- `turbine_ct::Array{TF,nTurbines}`: turbine thrust coefficients for the given state
- `turbine_ai::Array{TF,nTurbines}`: turbine axial induction for the given state
- `rotor_diameter::Array{TF,nTurbines}`: turbine rotor diameters
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_local_ti::Array{TF,nTurbines}`: turbine local turbulence intensity for 
    the given state
- `sorted_turbine_index::Array{TF,nTurbines}`: array containing indices of wind turbines 
    from most upwind to most downwind turbine in the given state
- `wtvelocities::Array{TF,nTurbines}`: effective inflow wind speed for given state
- `wind_resource::DiscretizedWindResource`: contains wind resource discreption (directions,
    speeds, frequencies, etc)
- `wind_farm_state_id::Int`: index to correct state to use from wind resource provided.
    Defaults to 1
- `downwind_turbine_id::Int`: index of wind turbine of interest (if any). If not a point for
    calculating effective wind speed of a turbine, then provide 0 (default)
"""
function point_velocity(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
                    wind_resource, model_set::AbstractModelSet;
                    wind_farm_state_id=1, downwind_turbine_id=0)

    wakedeficitmodel = model_set.wake_deficit_model
    wakedeflectionmodel = model_set.wake_deflection_model
    wakecombinationmodel = model_set.wake_combination_model
    
    # extract flow information
    wind_speed = wind_resource.wind_speeds[wind_farm_state_id]
    reference_height = wind_resource.measurement_heights[wind_farm_state_id]

    # set ground height 
    ground_height = wind_resource.wind_shear_model.ground_height    # TODO: allow topology to be given

    # find order for wind shear and deficit calculations
    shear_order = wind_resource.wind_shear_model.shear_order

    # adjust wind speed for wind shear
    if shear_order == "nothing"
        wind_speed_internal = wind_speed
    elseif shear_order == "first"
        wind_speed_internal = adjust_for_wind_shear(locz, wind_speed, reference_height, ground_height, wind_resource.wind_shear_model)
    else
        wind_speed_internal = wind_speed
    end

    # get number of turbines
    nturbines = length(turbine_x)

    # initialize deficit summation term to zero
    deficit_sum = 0.0

    # loop through all turbines
    for u=1:nturbines

        # get index of upstream turbine
        upwind_turb_id = Int(sorted_turbine_index[u])

        # don't allow turbine to impact itself
        if upwind_turb_id == downwind_turbine_id; continue; end

        # downstream distance between upstream turbine and point
        x = locx - turbine_x[upwind_turb_id]

        # check turbine relative locations
        if x > 1E-6
            # skip this loop if it would include a turbine's impact on itself)
            if upwind_turb_id==downwind_turbine_id; continue; end

            # calculate wake deflection of the current wake at the point of interest
            horizontal_deflection = wake_deflection_model(locx, locy, locz, turbine_x, turbine_yaw, turbine_ct,
                            upwind_turb_id, rotor_diameter, turbine_local_ti, wakedeflectionmodel)

            vertical_deflection = 0.0

            # velocity difference in the wake
            deltav = wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, horizontal_deflection, vertical_deflection,
                            upwind_turb_id, downwind_turbine_id, hub_height, rotor_diameter, turbine_ai,
                            turbine_local_ti, turbine_ct, turbine_yaw, wakedeficitmodel)

            # combine deficits according to selected wake combination method
            deficit_sum = wake_combination_model(deltav, wind_speed_internal, wtvelocities[upwind_turb_id], deficit_sum, wakecombinationmodel)
            # println(deficit_sum, " ", downwind_turbine_id, " ", upwind_turb_id)
            
        end
    end

    # find velocity at point without shear
    point_velocity = wind_speed_internal - deficit_sum

    if shear_order == "nothing"
        point_velocity_out = point_velocity
    elseif shear_order == "first"
        point_velocity_out = point_velocity
    else
        point_velocity_out = adjust_for_wind_shear(locz, point_velocity, reference_height, ground_height, wind_resource.wind_shear_model)        
    end

    return point_velocity_out

end

"""
    wake_added_tilt(loc, turbine_x, turbine_y, turbine_z, turbine_tilt, TSR, turbine_ct, turbine_ai,
    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
    wind_resource, model_set::AbstractModelSet;
    wind_farm_state_id=1, downwind_turbine_id=0)

Calculates the added tilt due to secondary wake steering

# Arguments
- `U_inf::Float`: streamwise wind velocity at turbine being compared to
- `W::Float`: vertical spanwise wind velocity at turbine being compared to
- `U_inf_initial::Float`: streamwise wind velocity of wind farm
- `deltay::Float`: distance in y-direction between turbine of interest and turbine being compared to
- `z_i::Float`: vertical location of center of turbine being compared to
- `rotor_diameter::Float`: rotor diameter of turbine being comapared to
- `hub_height::Float`: hub height of turbine being compared to
- `cT::Float`: coefficient of thrust for turbine being compared to
- `axial_induction::Float`: axial induction factor for turbine being compared to
"""
function wake_added_tilt(U_inf, W, U_inf_initial, deltay, z_i, rotor_diameter, hub_height, cT, TSR, axial_induction)

    # turbine parameters
    D = rotor_diameter
    HH = hub_height
    ai = axial_induction
    avgW = W

    # flow parameters
    Uinf = U_inf_initial

    # epsilon gain
    eps_gain = 0.2
    eps = eps_gain * D

    # find velocity at top and bottom of rotor swept area
    vel_top = ((HH+D/2)/HH)^0.12
    vel_bottom = ((HH-D/2)/HH)^0.12

    # find Gamma at the top and bottom of the rotor swept area
    # Gamma_top = gamma(D, vel_top, Uinf, cT)
    Gamma_top = (pi/8)*D*Uinf*cT*vel_top                # why do we use Uinf here, should we use U_inf?
    # Gamma_bot = -1*gamma(D, vel_bottom, Uinf, cT)
    Gamma_bot = (pi/8)*D*Uinf*cT*vel_bottom*-1.0

    # Use turbine average velocity to find Gamma due to wake rotation
    turbine_average_velocity = U_inf
    Gamma_wake_rotation = pi*D*(ai-ai^2)*turbine_average_velocity/TSR

    ### Calculate the spanwise and vertical velocities induced by tilt ###
    # top vortex
    dz_top = z_i - (HH+D/2)             
    r_top = deltay^2 + dz_top^2             
    core_shape = 1-exp(-r_top/(eps^2))          
    v_top = ((Gamma_top*dz_top)/(2*pi*r_top))*core_shape

    # bottom vortex
    dz_bot = z_i - (HH - D/2)
    r_bot = deltay^2 + dz_bot^2
    core_shape = 1-exp(-r_bot/(eps^2))
    v_bot = ((Gamma_bot*dz_bot)/(2*pi*r_bot))*core_shape

    # wake rotation vortex
    dz_center = z_i - HH
    r_center = deltay^2 + dz_center^2
    core_shape = 1 - exp(-r_center/(eps^2))
    v_center = ((Gamma_wake_rotation*dz_center)/(2*pi*r_center))*core_shape

    val = 2*(avgW-v_center)/(v_top+v_bot)           # why is this multiplied by 2?

    # cap the added_tilt to be between -45 and 45
    if val > 1.0
        val = 1.0
    elseif val < -1.0
        val = -1.0
    end

    added_tilt = 0.5*asin(val)              # why is this multiplied by 0.5?

    return added_tilt
end

"""
    point_velocity_tilt(loc, turbine_x, turbine_y, turbine_z, turbine_tilt, TSR, turbine_ct, turbine_ai,
    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
    wind_resource, model_set::AbstractModelSet;
    wind_farm_state_id=1, downwind_turbine_id=0)

Calculates the wind speed at a given point for a given state (with tilt)

# Arguments
- `loc::Array{TF,3}`: Location of interest
- `turbine_x::Array{TF,nTurbines}`: turbine east-west locations in the state 
    reference frame
- `turbine_y::Array{TF,nTurbines}`: turbine north-south locations in the state 
    reference frame
- `turbine_z::Array{TF,nTurbines}`: turbine base height in the state reference frame
- `turbine_tilt::Array{TF,nTurbines}`: turbine tilt for the given wind direction in 
    radians
- `turbine_ct::Array{TF,nTurbines}`: turbine thrust coefficients for the given state
- `turbine_ai::Array{TF,nTurbines}`: turbine axial induction for the given state
- `rotor_diameter::Array{TF,nTurbines}`: turbine rotor diameters
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_local_ti::Array{TF,nTurbines}`: turbine local turbulence intensity for 
    the given state
- `sorted_turbine_index::Array{TF,nTurbines}`: array containing indices of wind turbines 
    from most upwind to most downwind turbine in the given state
- `wtvelocities::Array{TF,nTurbines}`: effective inflow wind speed for given state
- `wind_resource::DiscretizedWindResource`: contains wind resource discreption (directions,
    speeds, frequencies, etc)
- `wind_farm_state_id::Int`: index to correct state to use from wind resource provided.
    Defaults to 1
- `downwind_turbine_id::Int`: index of wind turbine of interest (if any). If not a point for
    calculating effective wind speed of a turbine, then provide 0 (default)
"""
function point_velocity_tilt(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_tilt, TSR, turbine_ct, turbine_ai,
                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
                    wind_resource, model_set::AbstractModelSet;
                    wind_farm_state_id=1, downwind_turbine_id=0)

    wakedeficitmodel = model_set.wake_deficit_model
    wakedeflectionmodel = model_set.wake_deflection_model
    wakecombinationmodel = model_set.wake_combination_model
    
    # extract flow information
    wind_speed = wind_resource.wind_speeds[wind_farm_state_id]
    reference_height = wind_resource.measurement_heights[wind_farm_state_id]

    # set ground height 
    ground_height = wind_resource.wind_shear_model.ground_height    # TODO: allow topology to be given

    # find order for wind shear and deficit calculations
    shear_order = wind_resource.wind_shear_model.shear_order

    # adjust wind speed for wind shear
    if shear_order == "nothing"
        wind_speed_internal = wind_speed
    elseif shear_order == "first"
        wind_speed_internal = adjust_for_wind_shear(locz, wind_speed, reference_height, ground_height, wind_resource.wind_shear_model)
    else
        wind_speed_internal = wind_speed
    end

    # get number of turbines
    nturbines = length(turbine_x)

    # initialize deficit summation term to zero
    deficit_sum = 0.0

    # loop through all turbines
    for u=1:nturbines

        # get index of upstream turbine
        upwind_turb_id = Int(sorted_turbine_index[u])

        # don't allow turbine to impact itself
        if upwind_turb_id == downwind_turbine_id; continue; end

        # downstream distance between upstream turbine and point
        x = locx - turbine_x[upwind_turb_id]

        # check turbine relative locations
        if x > 1E-6
            # skip this loop if it would include a turbine's impact on itself)
            if upwind_turb_id==downwind_turbine_id; continue; end

            # find the added tilt angle due to the vortices
            added_tilt = wake_added_tilt()
            turbine_tilt += added_tilt

            # calculate wake deflection of the current wake at the point of interest
            horizontal_deflection = 0.0

            vertical_deflection = wake_deflection_model(locx, locy, locz, turbine_x, turbine_tilt, turbine_ct,
            upwind_turb_id, rotor_diameter, turbine_local_ti, wakedeflectionmodel)

            # velocity difference in the wake
            deltav = wake_deficit_model(locx, locy, locz, turbine_x, turbine_y, turbine_z, horizontal_deflection, vertical_deflection,
                            upwind_turb_id, downwind_turbine_id, hub_height, rotor_diameter, turbine_ai,
                            turbine_local_ti, turbine_ct, turbine_tilt, wakedeficitmodel)

            # combine deficits according to selected wake combination method
            deficit_sum = wake_combination_model(deltav, wind_speed_internal, wtvelocities[upwind_turb_id], deficit_sum, wakecombinationmodel)
            # println(deficit_sum, " ", downwind_turbine_id, " ", upwind_turb_id)
            
        end
    end

    # find velocity at point without shear
    point_velocity = wind_speed_internal - deficit_sum

    if shear_order == "nothing"
        point_velocity_out = point_velocity
    elseif shear_order == "first"
        point_velocity_out = point_velocity
    else
        point_velocity_out = adjust_for_wind_shear(locz, point_velocity, reference_height, ground_height, wind_resource.wind_shear_model)        
    end

    return point_velocity_out

end

"""
    point_velocity(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
    sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
    model_set::AbstractModelSet; wind_farm_state_id=1)

Calculates the wind speed at a given point for a given state

# Arguments
- `turbine_x::Array{TF,nTurbines}`: turbine east-west locations in the state 
    reference frame
- `turbine_y::Array{TF,nTurbines}`: turbine north-south locations in the state 
    reference frame
- `turbine_z::Array{TF,nTurbines}`: turbine base height in the state reference frame
- `rotor_diameter::Array{TF,nTurbines}`: turbine rotor diameters
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_yaw::Array{TF,nTurbines}`: turbine yaw for the given wind direction in 
    radians
- `sorted_turbine_index::Array{TF,nTurbines}`: turbine sorted order upstream to downstream 
    for given state
- `ct_model::AbstractThrustCoefficientModel`: defines how the thrust coefficient changes 
    with state etc
- rotor_sample_points_y::Array{TF,N}`: horizontal wind location of points to sample across 
    the rotor swept area when calculating the effective wind speed for the wind turbine. 
    Points are centered at the hub (0,0) and scaled by the radius (1=tip of blades) 
- rotor_sample_points_z::Array{TF,N}`: vertical wind location of points to sample across the 
    rotor swept area when calculating the effective wind speed for the wind turbine. Points
    are centered at the hub (0,0) and scaled by the radius (1=tip of blades)
- `wind_resource::DiscretizedWindResource`: wind resource discreption (directions, speeds, 
    frequencies, etc)
- `model_set::AbstractModelSet`: defines wake-realated models to be used in analysis
- `wind_farm_state_id::Int`: index to correct state to use from wind resource provided.
    Defaults to 1
"""
function turbine_velocities_one_direction(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
                    sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set::AbstractModelSet; wind_farm_state_id=1, velocity_only=true)
    
    # get number of turbines and rotor sample point
    n_turbines = length(turbine_x)
    n_rotor_sample_points = length(rotor_sample_points_y)

    # initialize correct array types
    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                            typeof(hub_height[1]),typeof(turbine_yaw[1]))

    # initialize arrays
    turbine_velocities = zeros(arr_type, n_turbines)
    turbine_ct = zeros(arr_type, n_turbines)
    turbine_ai = zeros(arr_type, n_turbines)
    turbine_local_ti = zeros(arr_type, n_turbines)

    # loop over all turbines
    for d=1:n_turbines

        # get index of downstream turbine
        downwind_turbine_id = Int(sorted_turbine_index[d])

        # initialize downstream wind turbine velocity to zero
        # println("start array: ", turbine_velocities[downwind_turbine_id])
        # wind_turbine_velocity = typeof(turbine_velocities[downwind_turbine_id])(0.0)
        wind_turbine_velocity = 0.0
        # turbine_velocities[downwind_turbine_id] = 0.0

        # loop over all rotor sample points to approximate the effective inflow velocity
        for p=1:n_rotor_sample_points

            # scale rotor sample point coordinate by rotor diameter (in rotor hub ref. frame)
            local_rotor_sample_point_y = rotor_sample_points_y[p]*0.5*rotor_diameter[downwind_turbine_id]
            local_rotor_sample_point_z = rotor_sample_points_z[p]*0.5*rotor_diameter[downwind_turbine_id]

            # put rotor sample points in wind direction coordinate system, and account for yaw
            locx = turbine_x[downwind_turbine_id] .+ local_rotor_sample_point_y*sin(turbine_yaw[downwind_turbine_id])
            locy = turbine_y[downwind_turbine_id] .+ local_rotor_sample_point_y*cos(turbine_yaw[downwind_turbine_id])
            locz = turbine_z[downwind_turbine_id] .+ hub_height[downwind_turbine_id] + local_rotor_sample_point_z

            # calculate the velocity at given point
            point_velocity_with_shear = point_velocity(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
                                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, turbine_velocities,
                                    wind_resource, model_set,
                                    wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=downwind_turbine_id)

            # add sample point velocity to turbine velocity to be averaged later
            wind_turbine_velocity += point_velocity_with_shear

        end

        # final velocity calculation for downstream turbine (average equally across all points)
        wind_turbine_velocity /= n_rotor_sample_points

        turbine_velocities[downwind_turbine_id] = deepcopy(wind_turbine_velocity)

        # update thrust coefficient for downstream turbine
        turbine_ct[downwind_turbine_id] = calculate_ct(turbine_velocities[downwind_turbine_id], ct_model[downwind_turbine_id])

        # update axial induction for downstream turbine
        turbine_ai[downwind_turbine_id] = _ct_to_axial_ind_func(turbine_ct[downwind_turbine_id], turbine_yaw[downwind_turbine_id])

        # get local turbulence intensity for this wind state
        ambient_ti = wind_resource.ambient_tis[wind_farm_state_id]
        
        # update local turbulence intensity for downstream turbine
        turbine_local_ti[downwind_turbine_id] = calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
                            turbine_velocities, turbine_ct, model_set.local_ti_model; turbine_id=downwind_turbine_id, tol=1E-6)

        # println("local ti turb 9: ", turbine_local_ti[downwind_turbine_id])
    end

    # df = DataFrame(ID=1:n_turbines, V=turbine_velocities, TI=turbine_local_ti, CT=turbine_ct)
    # CSV.write("internaldata.txt", df)

    if velocity_only
        return turbine_velocities 
    else
        return turbine_velocities, turbine_ct, turbine_ai, turbine_local_ti
    end
end

"""
    turbine_velocities_one_direction_tilt(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_tilt, TSR,
    sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
    model_set::AbstractModelSet; wind_farm_state_id=1)

Calculates the wind speed at a given point for a given state

# Arguments
- `turbine_x::Array{TF,nTurbines}`: turbine east-west locations in the state 
    reference frame
- `turbine_y::Array{TF,nTurbines}`: turbine north-south locations in the state 
    reference frame
- `turbine_z::Array{TF,nTurbines}`: turbine base height in the state reference frame
- `rotor_diameter::Array{TF,nTurbines}`: turbine rotor diameters
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_tilt::Array{TF,nTurbines}`: turbine tilt for the given wind direction in 
    radians
- `sorted_turbine_index::Array{TF,nTurbines}`: turbine sorted order upstream to downstream 
    for given state
- `ct_model::AbstractThrustCoefficientModel`: defines how the thrust coefficient changes 
    with state etc
- rotor_sample_points_y::Array{TF,N}`: horizontal wind location of points to sample across 
    the rotor swept area when calculating the effective wind speed for the wind turbine. 
    Points are centered at the hub (0,0) and scaled by the radius (1=tip of blades) 
- rotor_sample_points_z::Array{TF,N}`: vertical wind location of points to sample across the 
    rotor swept area when calculating the effective wind speed for the wind turbine. Points
    are centered at the hub (0,0) and scaled by the radius (1=tip of blades)
- `wind_resource::DiscretizedWindResource`: wind resource discreption (directions, speeds, 
    frequencies, etc)
- `model_set::AbstractModelSet`: defines wake-realated models to be used in analysis
- `wind_farm_state_id::Int`: index to correct state to use from wind resource provided.
    Defaults to 1
"""
function turbine_velocities_one_direction_tilt(turbine_x, turbine_y, turbine_z, rotor_diameter, hub_height, turbine_tilt, TSR,
                    sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                    model_set::AbstractModelSet; wind_farm_state_id=1, velocity_only=true)
    
    # get number of turbines and rotor sample point
    n_turbines = length(turbine_x)
    n_rotor_sample_points = length(rotor_sample_points_y)

    # initialize correct array types
    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                            typeof(hub_height[1]),typeof(turbine_tilt[1]))

    # initialize arrays
    turbine_velocities = zeros(arr_type, n_turbines)
    turbine_ct = zeros(arr_type, n_turbines)
    turbine_ai = zeros(arr_type, n_turbines)
    turbine_local_ti = zeros(arr_type, n_turbines)

    # loop over all turbines
    for d=1:n_turbines

        # get index of downstream turbine
        downwind_turbine_id = Int(sorted_turbine_index[d])

        # initialize downstream wind turbine velocity to zero
        # println("start array: ", turbine_velocities[downwind_turbine_id])
        # wind_turbine_velocity = typeof(turbine_velocities[downwind_turbine_id])(0.0)
        wind_turbine_velocity = 0.0
        # turbine_velocities[downwind_turbine_id] = 0.0

        # loop over all rotor sample points to approximate the effective inflow velocity
        for p=1:n_rotor_sample_points

            # scale rotor sample point coordinate by rotor diameter (in rotor hub ref. frame)
            local_rotor_sample_point_y = rotor_sample_points_y[p]*0.5*rotor_diameter[downwind_turbine_id]
            local_rotor_sample_point_z = rotor_sample_points_z[p]*0.5*rotor_diameter[downwind_turbine_id]

            # put rotor sample points in wind direction coordinate system, and account for yaw
            locx = turbine_x[downwind_turbine_id] .+ local_rotor_sample_point_y*sin(turbine_tilt[downwind_turbine_id])
            locy = turbine_y[downwind_turbine_id] .+ local_rotor_sample_point_y
            locz = turbine_z[downwind_turbine_id] .+ hub_height[downwind_turbine_id] .+ local_rotor_sample_point_z*cos(turbine_tilt[downwind_turbine_id])

            # calculate the velocity at given point
            point_velocity_with_shear = point_velocity_tilt(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_tilt, TSR, turbine_ct, turbine_ai,
                                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, turbine_velocities,
                                    wind_resource, model_set,
                                    wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=downwind_turbine_id)

            # add sample point velocity to turbine velocity to be averaged later
            wind_turbine_velocity += point_velocity_with_shear

        end

        # final velocity calculation for downstream turbine (average equally across all points)
        wind_turbine_velocity /= n_rotor_sample_points

        turbine_velocities[downwind_turbine_id] = deepcopy(wind_turbine_velocity)

        # update thrust coefficient for downstream turbine
        turbine_ct[downwind_turbine_id] = calculate_ct(turbine_velocities[downwind_turbine_id], ct_model[downwind_turbine_id])

        # update axial induction for downstream turbine
        turbine_ai[downwind_turbine_id] = _ct_to_axial_ind_func(turbine_ct[downwind_turbine_id], turbine_yaw[downwind_turbine_id])

        # get local turbulence intensity for this wind state
        ambient_ti = wind_resource.ambient_tis[wind_farm_state_id]
        
        # update local turbulence intensity for downstream turbine
        turbine_local_ti[downwind_turbine_id] = calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_tilt, turbine_local_ti, sorted_turbine_index,
                            turbine_velocities, turbine_ct, model_set.local_ti_model; turbine_id=downwind_turbine_id, tol=1E-6)

        # println("local ti turb 9: ", turbine_local_ti[downwind_turbine_id])
    end

    # df = DataFrame(ID=1:n_turbines, V=turbine_velocities, TI=turbine_local_ti, CT=turbine_ct)
    # CSV.write("internaldata.txt", df)

    if velocity_only
        return turbine_velocities 
    else
        return turbine_velocities, turbine_ct, turbine_ai, turbine_local_ti
    end
end

function turbine_velocities_one_direction(x, turbine_z, rotor_diameter, hub_height, turbine_yaw,
    sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
    model_set::AbstractModelSet; wind_farm_state_id=1, velocity_only=true)

    n_turbines = Int(length(x)/2)
    # println(typeof(x), n_turbines)
    turbine_x = x[1:n_turbines] 
    turbine_y = x[n_turbines+1:end]
    # println(turbine_x)
    # println("turbine_x type ", typeof(turbine_x))
    # println("type of x ", typeof(x))

    # get number of turbines and rotor sample point
    # n_turbines = length(turbine_x)
    n_rotor_sample_points = length(rotor_sample_points_y)

    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),
                typeof(hub_height[1]),typeof(turbine_yaw[1]))
    turbine_velocities = zeros(arr_type, n_turbines)
    turbine_ct = zeros(arr_type, n_turbines)
    turbine_ai = zeros(arr_type, n_turbines)
    turbine_local_ti = zeros(arr_type, n_turbines)

    for d=1:n_turbines

        # get index of downstream turbine
        downwind_turbine_id = Int(sorted_turbine_index[d])

        # initialize downstream wind turbine velocity to zero
        # println("start array: ", turbine_velocities[downwind_turbine_id])
        # wind_turbine_velocity = typeof(turbine_velocities[downwind_turbine_id])(0.0)
        wind_turbine_velocity = 0.0
        # turbine_velocities[downwind_turbine_id] = 0.0

        for p=1:n_rotor_sample_points


            # scale rotor sample point coordinate by rotor diameter (in rotor hub ref. frame)
            local_rotor_sample_point_y = rotor_sample_points_y[p]*0.5*rotor_diameter[downwind_turbine_id]
            local_rotor_sample_point_z = rotor_sample_points_z[p]*0.5*rotor_diameter[downwind_turbine_id]

            locx = turbine_x[downwind_turbine_id] .+ local_rotor_sample_point_y*sin(turbine_yaw[downwind_turbine_id])
            locy = turbine_y[downwind_turbine_id] .+ local_rotor_sample_point_y*cos(turbine_yaw[downwind_turbine_id])
            locz = turbine_z[downwind_turbine_id] .+ hub_height[downwind_turbine_id] + local_rotor_sample_point_z

            # calculate the velocity at given point
            point_velocity_with_shear = point_velocity(locx, locy, locz, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
                                rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, turbine_velocities,
                                wind_resource, model_set,
                                wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=downwind_turbine_id)

            # add sample point velocity to turbine velocity to be averaged later
            wind_turbine_velocity += point_velocity_with_shear

        end

        # final velocity calculation for downstream turbine (average equally across all points)
        wind_turbine_velocity /= n_rotor_sample_points

        turbine_velocities[downwind_turbine_id] = deepcopy(wind_turbine_velocity)

        # update thrust coefficient for downstream turbine
        turbine_ct[downwind_turbine_id] = calculate_ct(turbine_velocities[downwind_turbine_id], ct_model[downwind_turbine_id])

        # update axial induction for downstream turbine
        turbine_ai[downwind_turbine_id] = _ct_to_axial_ind_func(turbine_ct[downwind_turbine_id], turbine_yaw[downwind_turbine_id])

        # update local turbulence intensity for downstream turbine
        ambient_ti = wind_resource.ambient_tis[wind_farm_state_id]
        turbine_local_ti[downwind_turbine_id] = calculate_local_ti(turbine_x, turbine_y, ambient_ti, rotor_diameter, hub_height, turbine_yaw, turbine_local_ti, sorted_turbine_index,
                    turbine_velocities, turbine_ct, model_set.local_ti_model; turbine_id=downwind_turbine_id, tol=1E-6)

    end

    if velocity_only
        return turbine_velocities 
    else
        return turbine_velocities, turbine_ct, turbine_ai, turbine_local_ti
    end

end

# turbine_velocities_one_direction!(model_set::AbstractModelSet, problem_description::AbstractWindFarmProblem; wind_farm_state_id=1) = turbine_velocities_one_direction!([0.0], [0.0],
# model_set::AbstractModelSet, problem_description::AbstractWindFarmProblem; wind_farm_state_id=1)

"""
calculate_flow_field(xrange, yrange, zrange, model_set::AbstractModelSet, turbine_x, 
    turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai, rotor_diameter, hub_height, 
    turbine_local_ti, sorted_turbine_index, wtvelocities, wind_resource; wind_farm_state_id=1)

Generates a flow field for a given state and cross section

# Arguments
- `xrange::Range`: range defining east-west locations to sample in global reference frame
- `yrange::Range`: range defining north-west locations to sample in global reference frame
- `zrange::Range`: range defining vertical locations to sample in global reference frame
- `model_set::AbstractModelSet`: defines wake-realated models to be used in analysis
- `turbine_x::Array{TF,nTurbines}`: turbine east-west locations in the global 
    reference frame
- `turbine_y::Array{TF,nTurbines}`: turbine north-south locations in the global 
    reference frame
- `turbine_z::Array{TF,nTurbines}`: turbine base height in the global reference frame
- `turbine_yaw::Array{TF,nTurbines}`: turbine yaw for the given wind direction in 
    radians
- `turbine_ct::Array{TF,nTurbines}`: thrust coefficient of each turbine for the given state
- `turbine_ai::Array{TF,nTurbines}`: turbine axial induction for the given state
- `rotor_diameter::Array{TF,nTurbines}`: turbine rotor diameters
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_local_ti::Array{TF,nTurbines}`: turbine local turbulence intensity for 
    the given state
- `sorted_turbine_index::Array{TF,nTurbines}`: turbine north-south locations in the 
    global reference frame
- `wtvelocities::Array{TF,nTurbines}`: effective inflow wind speed for given state
- `wind_resource::DiscretizedWindResource`: wind resource discreption (directions, speeds, 
    frequencies, etc)
- `wind_farm_state_id::Int`: index to correct state to use from wind resource provided.
    Defaults to 1
"""
function calculate_flow_field(xrange, yrange, zrange,
    model_set::AbstractModelSet, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
    wind_resource; wind_farm_state_id=1)

    xlen = length(xrange)
    ylen = length(yrange)
    zlen = length(zrange)
    npoints = xlen*ylen*zlen
    point_velocities = zeros(npoints)
    point_velocities = reshape(point_velocities, (zlen, ylen, xlen))

    # rotate to direction frame for velocity calculations
    rot_tx, rot_ty = rotate_to_wind_direction(turbine_x, turbine_y, wind_resource.wind_directions[wind_farm_state_id])

    # sort the turbines
    sorted_turbine_index = sortperm(rot_tx)

    for zi in 1:zlen
        for yi in 1:ylen
            for xi in 1:xlen
                locx = xrange[xi]
                locy = yrange[yi]
                locz = zrange[zi]
                locx, locy = rotate_to_wind_direction(locx, locy, wind_resource.wind_directions[wind_farm_state_id])

                point_velocities[zi, yi, xi] = point_velocity(locx, locy, locz, rot_tx, rot_ty, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
                    wind_resource, model_set,
                    wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=0)

            end
        end
    end

    

    # if zlen == 1
    #     return point_velocities[1,1:ylen,1:xlen]
    # elseif ylen == 1
    #     return point_velocities[1:zlen,1,1:xlen]
    # elseif xlen == 1
    #     return point_velocities[1:zlen,1:ylen,1]
    # else
    return point_velocities[1:zlen,1:ylen,1:xlen]
    # end

end

function calculate_flow_field(xrange, yrange, zrange,
    model_set::AbstractModelSet, turbine_x, turbine_y, turbine_z, turbine_yaw,
    rotor_diameter, hub_height, ct_models, rotor_sample_points_y, rotor_sample_points_z,
    wind_resource; wind_farm_state_id=1)

    # rotate to direction frame for velocity calculations
    rot_tx, rot_ty = rotate_to_wind_direction(turbine_x, turbine_y, wind_resource.wind_directions[wind_farm_state_id])

    # sort the turbines
    sorted_turbine_index = sortperm(rot_tx)

    turbine_velocities, turbine_ct, turbine_ai, turbine_local_ti = turbine_velocities_one_direction(rot_tx, rot_ty, turbine_z, rotor_diameter, hub_height, turbine_yaw,
    sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
    model_set, wind_farm_state_id=wind_farm_state_id, velocity_only=false)

    return calculate_flow_field(xrange, yrange, zrange,
        model_set, turbine_x, turbine_y, turbine_z, turbine_yaw, turbine_ct, turbine_ai,
        rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, turbine_velocities,
        wind_resource, wind_farm_state_id=wind_farm_state_id)

end

"""
calculate_flow_field_tilt(xrange, yrange, zrange, model_set::AbstractModelSet, turbine_x, 
    turbine_y, turbine_z, turbine_tilt, turbine_ct, turbine_ai, rotor_diameter, hub_height, 
    turbine_local_ti, sorted_turbine_index, wtvelocities, wind_resource; wind_farm_state_id=1)

Generates a flow field for a given state and cross section

# Arguments
- `xrange::Range`: range defining east-west locations to sample in global reference frame
- `yrange::Range`: range defining north-west locations to sample in global reference frame
- `zrange::Range`: range defining vertical locations to sample in global reference frame
- `model_set::AbstractModelSet`: defines wake-realated models to be used in analysis
- `turbine_x::Array{TF,nTurbines}`: turbine east-west locations in the global 
    reference frame
- `turbine_y::Array{TF,nTurbines}`: turbine north-south locations in the global 
    reference frame
- `turbine_z::Array{TF,nTurbines}`: turbine base height in the global reference frame
- `turbine_tilt::Array{TF,nTurbines}`: turbine yaw for the given wind direction in 
    radians
- `turbine_ct::Array{TF,nTurbines}`: thrust coefficient of each turbine for the given state
- `turbine_ai::Array{TF,nTurbines}`: turbine axial induction for the given state
- `rotor_diameter::Array{TF,nTurbines}`: turbine rotor diameters
- `hub_height::Array{TF,nTurbines}`: turbine hub heights
- `turbine_local_ti::Array{TF,nTurbines}`: turbine local turbulence intensity for 
    the given state
- `sorted_turbine_index::Array{TF,nTurbines}`: turbine north-south locations in the 
    global reference frame
- `wtvelocities::Array{TF,nTurbines}`: effective inflow wind speed for given state
- `wind_resource::DiscretizedWindResource`: wind resource discreption (directions, speeds, 
    frequencies, etc)
- `wind_farm_state_id::Int`: index to correct state to use from wind resource provided.
    Defaults to 1
"""

function calculate_flow_field_tilt(xrange, yrange, zrange,
    model_set::AbstractModelSet, turbine_x, turbine_y, turbine_z, turbine_tilt, TSR, turbine_ct, turbine_ai,
    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
    wind_resource; wind_farm_state_id=1)

    xlen = length(xrange)
    ylen = length(yrange)
    zlen = length(zrange)
    npoints = xlen*ylen*zlen
    point_velocities = zeros(npoints)
    point_velocities = reshape(point_velocities, (zlen, ylen, xlen))

    # rotate to direction frame for velocity calculations
    rot_tx, rot_ty = rotate_to_wind_direction(turbine_x, turbine_y, wind_resource.wind_directions[wind_farm_state_id])

    # sort the turbines
    sorted_turbine_index = sortperm(rot_tx)
    
    for zi in 1:zlen
        for yi in 1:ylen
            for xi in 1:xlen
                locx = xrange[xi]
                locy = yrange[yi]
                locz = zrange[zi]
                locx, locy = rotate_to_wind_direction(locx, locy, wind_resource.wind_directions[wind_farm_state_id])

                point_velocities[zi, yi, xi] = point_velocity_tilt(locx, locy, locz, rot_tx, rot_ty, turbine_z, turbine_tilt, TSR, turbine_ct, turbine_ai,
                    rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, wtvelocities,
                    wind_resource, model_set,
                    wind_farm_state_id=wind_farm_state_id, downwind_turbine_id=0)

            end
        end
    end

    

    # if zlen == 1
    #     return point_velocities[1,1:ylen,1:xlen]
    # elseif ylen == 1
    #     return point_velocities[1:zlen,1,1:xlen]
    # elseif xlen == 1
    #     return point_velocities[1:zlen,1:ylen,1]
    # else
    return point_velocities[1:zlen,1:ylen,1:xlen]
    # end

end

function calculate_flow_field_tilt(xrange, yrange, zrange,
    model_set::AbstractModelSet, turbine_x, turbine_y, turbine_z, turbine_tilt,
    rotor_diameter, hub_height, ct_models, TSR, rotor_sample_points_y, rotor_sample_points_z,
    wind_resource; wind_farm_state_id=1)

    # rotate to direction frame for velocity calculations
    rot_tx, rot_ty = rotate_to_wind_direction(turbine_x, turbine_y, wind_resource.wind_directions[wind_farm_state_id])

    # sort the turbines
    sorted_turbine_index = sortperm(rot_tx)

    turbine_velocities, turbine_ct, turbine_ai, turbine_local_ti = turbine_velocities_one_direction_tilt(rot_tx, rot_ty, turbine_z, rotor_diameter, hub_height, turbine_tilt, TSR,
    sorted_turbine_index, ct_models, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
    model_set, wind_farm_state_id=wind_farm_state_id, velocity_only=false)

    return calculate_flow_field_tilt(xrange, yrange, zrange,
        model_set, turbine_x, turbine_y, turbine_z, turbine_tilt, TSR, turbine_ct, turbine_ai,
        rotor_diameter, hub_height, turbine_local_ti, sorted_turbine_index, turbine_velocities,
        wind_resource, wind_farm_state_id=wind_farm_state_id)

end

