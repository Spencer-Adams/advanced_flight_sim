module sim_m
    use adams_m
    use jsonx_m
    use linalg_mod
    use micro_time_m
    use connection_m
    use vehicle_m
    implicit none
    ! Variables within sim_m
    type(vehicle_t), allocatable :: vehicles(:)
    integer :: num_vehicles 
    real :: init_state(13)
    real :: controls(4)
    real :: controls_from_connect(4)
    logical :: save_states, rk4_verbose
    type(connection) :: graphics, connect_controls 
    type(json_value), pointer :: j_main
    ! aero coefficients 
    integer :: newton_max_iter
contains
    ! end subroutine simulation_main
    subroutine run()
        implicit none 
        real :: y(13), y1(13), s(14)
        real :: dt, tf, t
        real :: cpu_start_time, cpu_end_time, time1, time2, actual_time, integrated_time
        integer :: io_unit, i 
        logical :: real_time
        ! initial conditions 
        t = 0.0
        y = 0.0
        real_time = .false.
        call jsonx_get(j_main, 'simulation.time_step[s]', dt, 0.0)
        call jsonx_get(j_main, 'simulation.end_time[s]', tf, 0.0)

        if (abs(dt)<TOLERANCE) then 
            real_time = .true.
            time1 = get_time()
            do i=1,num_vehicles 
                if(vehicles(i)%run_physics) call vehicle_tick_state(vehicles(i),t,dt)
            end do 
            time2 = get_time()

            dt = time2-time1 
            if (dt == 0.0) then 
                write(*,*) "THE DT IS EXACTLY ZERO!!!"
            end if 
            t = 0.0
        end if
        cpu_start_time = get_time()
        time1 = cpu_start_time
        integrated_time = 0.0
        do while(t<tf) ! while altitude is greater than 0 ft (altitude is positive going down in our coordinate systems) or time is less than final time
            ! write(*,*) "dt", dt
            ! if (is_use_controls) then 
                ! controls = connect_controls%recv()
            ! end if
            do i=1,num_vehicles 
                if(vehicles(i)%run_physics) call vehicle_tick_state(vehicles(i),t,dt)
            end do 
            t = t + dt 
            integrated_time = integrated_time + dt
            ! s(1) = t
            ! s(2:14) = y(1:13)
            ! write(*,'(14E22.13)') s
            ! call graphics%send(s)
            if(real_time) then ! get estimate for next dt time step
                time2 = get_time()
                dt = time2-time1 
                time1 = time2 
            end if 

        end do
        cpu_end_time = get_time()
        actual_time = cpu_end_time-cpu_start_time
        write(*,*) '    Total integrated time [s] = ', integrated_time
        write(*,*) 'Total actual elapsed time [s]= ', actual_time
        write(*,*) 'Total error in time [s] = ', integrated_time - actual_time
    end subroutine run

    subroutine init(filename)
        implicit none 
        character(100), intent(in) :: filename
        type(json_value), pointer :: j_connections, j_graphics, j_controls
        type(json_value), pointer :: j_vehicles, j_temp, j_atmosphere
        integer :: i 
        write(*,*) 'Initializing Simulation...'
        ! type2, intent(out) ::  arg2
        call jsonx_load(filename, j_main)
        ! global settings 
        ! write(*,*) 'Reading atmosphere object in json'
        ! call jsonx_get(j_main, 'atmosphere', j_atmosphere)

        call jsonx_get(j_main, 'simulation.rk4_verbose', rk4_verbose, .false.)
        call jsonx_get(j_main, 'simulation.save_states', save_states, .true.)
        call jsonx_get(j_main, 'simulation.geographic_model', geographic_model, 'none')
        geographic_model_ID = 0
        if (geographic_model == 'sphere') geographic_model_ID = 1
        if (geographic_model == 'ellipse') geographic_model_ID = 2
        
        write(*,*) 'Initializing vehicles'
        call jsonx_get(j_main, 'vehicles', j_vehicles)
        num_vehicles = json_value_count(j_vehicles)
        allocate(vehicles(num_vehicles))

        do i = 1, num_vehicles
            call json_value_get(j_vehicles, i, j_temp)
            ! call atmosphere_init(vehicles(i)%atm,j_atmosphere)
            call vehicle_init(vehicles(i), j_temp, save_states, rk4_verbose)
        end do 
    end subroutine init
end module sim_m
