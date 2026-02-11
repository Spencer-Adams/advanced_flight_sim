module vehicle_m
    use adams_m
    use jsonx_m
    use linalg_mod
    use micro_time_m
    use connection_m
    implicit none 

    character(len=:), allocatable :: geographic_model
    integer :: geographic_model_ID

    type vehicle_t
        type(json_value), pointer :: j_vehicle
        
        character(len=:), allocatable :: name
        character(len=:), allocatable :: type
        character(100) :: states_filename, rk4_filename, geographic_filename

        logical :: run_physics
        logical :: save_states, rk4_verbose
        integer :: iunit_states, iunit_rk4, iunit_trim, iunit_geographic
        real :: rho0

        ! mass constants
        real :: mass
        real :: I(3,3)
        real :: Ixxb, Iyyb, Izzb, Ixyb, Ixzb, Iyzb
        real :: Iinv(3,3)
        real :: h_gyro(3,3)
        real :: hdot_gyro(3)
        real :: hx, hy, hz
        real, allocatable :: h(:)

        ! aerodynamic constants
        real,allocatable :: aero_ref_location(:) ! has to be allocatable because will be read from json object
        real :: FM(6)
        real :: sref, long_ref, lat_ref
        real :: CL0, CLa, CLahat, CLqbar, CLde
        real :: CDL0, CDL1, CDL2, CDS2, CDqbar, CDaqbar, CDde, CDade, CDde2
        real :: CSb, CSpbar, CSapbar, CSrbar, CSda, CSdr
        real :: Cll0, Clb, Clpbar, Clrbar, Clarbar, Clda, Cldr
        real :: Cm0, Cma, Cmqbar, Cmahat, Cmde
        real :: Cnb, Cnpbar, Cnapbar, Cnrbar, Cnda, Cnada, Cndr
        real :: Thrust0, Ta, thrust_quat(4)
        real, allocatable :: thrust_location(:)

        ! stall model constants
        logical :: include_stall
        ! type(stall_settings_t) :: CLstall, CDstall, Cmstall
        real :: lambda_CL, lambda_CD, lambda_Cm
        real :: alpha_0CL, alpha_sCL
        real :: alpha_0CD, alpha_sCD
        real :: alpha_0Cm, alpha_sCm, Cmmin

        ! initialization constants
        real :: init_V, init_alt, init_state(13)
        real, allocatable :: init_eul(:) ! has to be allocatable because will be read from json object
        real :: latitude, longitude 

        ! variables
        real :: state(13)
        real :: controls(4)

        ! type(trim_settings_t) :: trim

    end type vehicle_t
contains 

    subroutine vehicle_init(this, j_vehicle_input, is_save_states, is_rk4_verbose)
        implicit none 
        type(vehicle_t) :: this 
        type(json_value), pointer :: j_vehicle_input
        logical, intent(in) :: is_save_states, is_rk4_verbose
        real :: denom 
        logical :: is_straight_fletchings
        character(len=:), allocatable :: init_type 
        real, allocatable :: thrust_orientation(:)
        real :: Z_temp,T_temp,P_temp,a_temp,mu_temp
        real :: euler_angles_init(3), azimuth_init
        this%j_vehicle => j_vehicle_input
        this%name = this%j_vehicle%name 
        write(*,*) ' Initializing ', this%name 
        call jsonx_get(this%j_vehicle, "type", this%type)
        write(*,*) '   - type = ', this%type 
        call jsonx_get(this%j_vehicle, 'run_physics', this%run_physics)
        this%save_states = is_save_states
        this%rk4_verbose = is_rk4_verbose
        ! get atmosphere stuff 
        call std_atm_English(0.0,Z_temp,T_temp,P_temp,this%rho0,a_temp,mu_temp)

        if(this%run_physics) then 
            if (this%save_states) then 
                this%states_filename = trim(this%name)//'_states.csv'
                open(newunit=this%iunit_states, file = this%states_filename, status = 'REPLACE')
                write(this%iunit_states,*) "t[s],u[ft/s],v[ft/s],"// &
                        "w[ft/s],p[rad/s],q[rad/s],"  // &
                        "r[rad/s],x[ft],y[ft]," // &                   
                        "z[ft],e0,ex," // &
                        "ey,ez"
                write(*,*) '   - saving states to ', this%states_filename
                
                if (geographic_model_ID > 0) then
                    this%geographic_filename = trim(this%name)//'_geographic.csv'
                    open(newunit=this%iunit_geographic, file = this%geographic_filename, status = 'REPLACE')
                    write(this%iunit_geographic,*) "time[s],latitude[deg],longitude[deg],azimuth[deg]"
                    write(*,*) '   - saving geographic data to ', this%geographic_filename
                end if
            end if 

            write(*,*) '   - Initializing mass and inertia properties'
            call mass_inertia(this)
            write(*,*) '   - Initializing aerodynamics'
            call jsonx_get(this%j_vehicle, "aerodynamics.reference.area[ft^2]", this%sref)
            call jsonx_get(this%j_vehicle, "aerodynamics.reference.longitudinal_length[ft]", this%long_ref)
            call jsonx_get(this%j_vehicle, "aerodynamics.reference.lateral_length[ft]", this%lat_ref)
            call jsonx_get(this%j_vehicle, "aerodynamics.reference.location[ft]", this%aero_ref_location,0.0,3)
            
            call jsonx_get(this%j_vehicle, "thrust.T0[lbf]", this%Thrust0,0.0)
            call jsonx_get(this%j_vehicle, "thrust.Ta", this%Ta,0.0)
            call jsonx_get(this%j_vehicle, "thrust.location[ft]", this%thrust_location,0.0,3)
            call jsonx_get(this%j_vehicle, "thrust.orientation[deg]", thrust_orientation, 0.0,3)
            thrust_orientation = thrust_orientation*PI/180.0
            this%thrust_quat = euler_to_quat(thrust_orientation)
            if (this%type == 'arrow') then
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CL.alpha", this%CLa)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.L0", this%CDL0)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.CL1_CL1", this%CDL2)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.is_straight_fletchings", is_straight_fletchings)
                if (is_straight_fletchings) then 
                    call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.0", this%Cll0)
                else
                    call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.1", this%Cll0)
                end if 
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.pbar", this%Clpbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.alpha", this%Cma)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.qbar", this%Cmqbar)
            end if
            if (this%type == 'aircraft') then 
                ! vehicle
                ! thrust
                ! coefficients 
                !CL 
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CL.0", this%CL0)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CL.alpha", this%CLa)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CL.alphahat", this%CLahat)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CL.qbar", this%CLqbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CL.elevator", this%CLde)
                !CS 
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.beta", this%CSb)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.pbar", this%CSpbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.alpha_pbar", this%CSapbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.rbar", this%CSrbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.aileron", this%CSda)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CS.rudder", this%CSdr)
                ! CD
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.L0", this%CDL0)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.CL1", this%CDL1)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.CL1_CL1", this%CDL2)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.CS_CS", this%CDS2)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.qbar", this%CDqbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.alpha_qbar", this%CDaqbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.elevator", this%CDde)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.alpha_elevator", this%CDade)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.elevator_elevator", this%CDde2)
                ! Cl
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.beta", this%Clb)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.pbar", this%Clpbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.rbar", this%Clrbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.alpha_rbar", this%Clarbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.aileron", this%Clda)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.rudder", this%Cldr)
                ! Cm 
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.0", this%Cm0)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.alpha", this%Cma)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.qbar", this%Cmqbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.alphahat", this%Cmahat)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cm.elevator", this%Cmde)
                ! Cn
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.beta", this%Cnb)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.pbar", this%Cnpbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.alpha_pbar", this%Cnapbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.rbar", this%Cnrbar)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.aileron", this%Cnda)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.alpha_aileron", this%Cnada)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cn.rudder", this%Cndr)
            end if 
            
            if(this%type == 'arrow' .or. this%type == 'aircraft') then 
                call jsonx_get(this%j_vehicle, 'aerodynamics.stall.include_stall', this%include_stall)
                if(this%include_stall) then
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.CL.alpha_0[deg]", this%alpha_0CL)
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.CL.alpha_s[deg]", this%alpha_sCL)
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.CL.lambda_b", this%lambda_CL)
                    this%alpha_0CL = this%alpha_0CL*PI/180.0
                    this%alpha_sCL = this%alpha_sCL*PI/180.0
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.CD.alpha_0[deg]", this%alpha_0CD)
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.CD.alpha_s[deg]", this%alpha_sCD)
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.CD.lambda_b", this%lambda_CD)
                    this%alpha_0CD = this%alpha_0CD*PI/180.0
                    this%alpha_sCD = this%alpha_sCD*PI/180.0
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.Cm.min", this%Cmmin)
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.Cm.alpha_0[deg]", this%alpha_0Cm)
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.Cm.alpha_s[deg]", this%alpha_sCm)
                    call jsonx_get(this%j_vehicle, "aerodynamics.stall.Cm.lambda_b", this%lambda_Cm)
                    this%alpha_0Cm = this%alpha_0Cm*PI/180.0
                    this%alpha_sCm = this%alpha_sCm*PI/180.0
                end if 
            end if
            
            this%init_state = 0.0
            ! call jsonx_get(this%j_vehicle, "initial.time[sec]", t, 0.0)
            call jsonx_get(this%j_vehicle, "initial.airspeed[ft/s]", this%init_V)
            write(*,*) "Initial velocity", this%init_V
            call jsonx_get(this%j_vehicle, "initial.altitude[ft]", this%init_state(9))
            this%init_state(9) = - this%init_state(9)
            this%init_alt = this%init_state(9)
            call jsonx_get(this%j_vehicle, "initial.latitude[deg]", this%latitude, 0.0)
            call jsonx_get(this%j_vehicle, "initial.longitude[deg]", this%longitude, 0.0)
            this%latitude = this%latitude*PI/180.0
            this%longitude = this%longitude*PI/180.0
            call jsonx_get(this%j_vehicle, "initial.Euler_angles[deg]", this%init_eul,0.0,3)
            this%init_eul = this%init_eul*PI/180.0
            call jsonx_get(this%j_vehicle, "initial.type", init_type)

            if (init_type=="state") then
                call init_to_state(this)
                this%init_state(10:13) = euler_to_quat(this%init_eul)
            else 
                call init_to_trim(this)
            end if 
            this%state = this%init_state

            call vehicle_write_state(this, 0.0, this%state)
            
            ! Write initial geographic data
            if (geographic_model_ID > 0) then
                euler_angles_init = quat_to_euler(this%state(10:13))
                azimuth_init = euler_angles_init(3)
                write(this%iunit_geographic,'(*(G0.15,:,","))') 0.0, this%latitude*180.0/PI, &
                      this%longitude*180.0/PI, azimuth_init*180.0/PI
            end if

        end if     

    end subroutine vehicle_init

    subroutine init_to_state(this)
        implicit none 
        type(vehicle_t) :: this
        type(json_value), pointer :: j_initial, j_state 
        real :: alpha, beta
        write(*,*) ' Setting State'
        call jsonx_get(this%j_vehicle, 'initial', j_initial)
        call jsonx_get(j_initial, 'state', j_state)
        call jsonx_get(j_state, 'angle_of_attack[deg]', alpha)
        alpha = alpha*PI/180.0 
        call jsonx_get(j_state, 'sideslip_angle[deg]', beta)
        beta = beta*PI/180.0
        this%init_state(1) = this%init_V*cos(alpha)*cos(beta)
        this%init_state(2) = this%init_V*sin(beta)
        this%init_state(3) = this%init_V*sin(alpha)*cos(beta)
        
        call jsonx_get(j_state, 'p[deg/s]', this%init_state(4))
        call jsonx_get(j_state, 'q[deg/s]', this%init_state(5))
        call jsonx_get(j_state, 'r[deg/s]', this%init_state(6))
        this%init_state(4:6) = this%init_state(4:6)*PI/180.0

        this%controls(:) = 0.0
        if (this%type == 'aircraft') then 
            call jsonx_get(j_state, 'aileron[deg]', this%controls(1))
            call jsonx_get(j_state, 'elevator[deg]', this%controls(2))
            call jsonx_get(j_state, 'rudder[deg]', this%controls(3))
            call jsonx_get(j_state, 'throttle', this%controls(4))
            this%controls(1:3) = this%controls(1:3)*PI/180.0
        end if 
    end subroutine init_to_state

    subroutine init_to_trim(this)
        implicit none 
        type(vehicle_t) :: this
        type(json_value), pointer :: j_initial, j_trim
        real :: alpha, beta, p_wind
        real :: trim_array(9)
        real :: trim_bank_angle, trim_elevation_angle, trim_sideslip_angle, trim_azimuth_angle
        real :: trim_load_factor, load_factor_temp, gravity_temp, a_c_temp
        real :: qt_temp(4), euler_temp(3), xdot_temp(3), v_t_temp
        real :: trim_climb_angle, finite_diff_step, relax_factor, newton_tol
        integer :: newton_max_iter
        logical :: is_trim_sideslip_angle, is_trim_sct_load_factor
        character(:), allocatable :: is_elevation_or_climb, trim_type, is_bank_or_beta_for_shss
        character(:), allocatable :: is_bank_or_load_factor_for_sct

        write(*,*) ' Setting trim'
        call jsonx_get(this%j_vehicle, 'initial', j_initial)
        call jsonx_get(j_initial, 'trim', j_trim)

        call jsonx_get(j_trim, "type", trim_type)
        call jsonx_get(j_trim, "elevation_or_climb", is_elevation_or_climb)
        if (is_elevation_or_climb == "climb") then
            call jsonx_get(j_trim, "climb_angle[deg]", trim_climb_angle)
            trim_climb_angle = trim_climb_angle*PI/180.0
            trim_elevation_angle = 0.0
        else 
            call jsonx_get(j_trim, "elevation_angle[deg]", trim_elevation_angle)
            trim_elevation_angle = trim_elevation_angle*PI/180.0
            trim_climb_angle = 0.0
        end if 
        call jsonx_get(j_trim, "solver.finite_difference_step_size", finite_diff_step)
        call jsonx_get(j_trim, "solver.relaxation_factor", relax_factor)
        call jsonx_get(j_trim, "solver.tolerance", newton_tol)
        call jsonx_get(j_trim, "solver.max_iterations", newton_max_iter)
        trim_bank_angle = 0.0
        trim_sideslip_angle = 0.0
        p_wind = 0.0
        trim_azimuth_angle = this%init_eul(3)
        is_trim_sideslip_angle = .false.
        is_trim_sct_load_factor = .false.
        trim_load_factor = 0.0
        if (trim_type == "sct") then
            call jsonx_get(j_trim, "type_sct.bank_or_load_factor", is_bank_or_load_factor_for_sct)
            if (is_bank_or_load_factor_for_sct == "bank") then 
                call jsonx_get(j_trim, "type_sct.bank_angle[deg]", trim_bank_angle)
            else if (is_bank_or_load_factor_for_sct == "load_factor") then 
                is_trim_sct_load_factor = .true.
                call jsonx_get(j_trim, "type_sct.load_factor", trim_load_factor)
            end if 
        else if (trim_type == "shss") then
            call jsonx_get(j_trim, "type_shss.bank_or_beta", is_bank_or_beta_for_shss)
            if (is_bank_or_beta_for_shss == "bank") then 
                call jsonx_get(j_trim, "type_shss.bank_angle[deg]", trim_bank_angle)
            else if (is_bank_or_beta_for_shss == "beta") then 
                is_trim_sideslip_angle = .true.
                call jsonx_get(j_trim, "type_shss.sideslip_angle[deg]", trim_sideslip_angle)
            end if 
        else if (trim_type == "vbr") then 
            call jsonx_get(j_trim, "type_vbr.p_wind[deg/s]", p_wind)
            call jsonx_get(j_trim, "type_vbr.bank_angle[deg]", trim_bank_angle)
        end if 
        trim_bank_angle = trim_bank_angle*PI/180.0
        trim_sideslip_angle = trim_sideslip_angle*PI/180.0
        p_wind = p_wind*PI/180.0
        alpha = 0.0
        beta = 0.0
        if (trim_type == "shss" .and. is_trim_sideslip_angle) then
            beta = trim_sideslip_angle
        end if 
        write(*,*) 
        write(*,*) 
        write(*,*) "alpha", alpha
        write(*,*) "beta", beta
        write(*,*) "x[ft]", this%init_state(7)
        write(*,*) "y[ft]", this%init_state(8)
        write(*,*) "altitude[ft]", this%init_state(9)
        write(*,*) "trim_type ", trim_type
        if (trim_type == "shss") then 
            write(*,*) "is_bank_or_beta_for_shss ", is_bank_or_beta_for_shss
        end if
        write(*,*) "trim_bank_angle", trim_bank_angle
        write(*,*) "trim_elevation_angle", trim_elevation_angle
        write(*,*) "trim_azimuth_angle", trim_azimuth_angle
        write(*,*) 
        write(*,*) 

        write(*,'(a)') 'Trimming Aircraft for '// trim_type
        write(*,'(a,f12.6)') '  --> Azimuth angle set to psi [deg] = ', trim_azimuth_angle*180.0/PI
        write(*,'(a,f12.6)') '  --> Elevation angle set to theta [deg] = ', trim_elevation_angle*180.0/PI
        write(*,'(a,f12.6)') '  --> Bank angle set to phi [deg] = ', trim_bank_angle*180.0/PI

        write(*,'(a,1x,e25.16)') 'Initial theta [deg] = ', trim_elevation_angle*180.0/PI
        write(*,'(a,1x,e25.16)') 'Initial gamma [deg] = ', trim_climb_angle * 180.0/PI
        write(*,'(a,1x,e25.16)') 'Initial phi [deg]   = ', trim_bank_angle*180.0/PI
        write(*,'(a,1x,e25.16)') 'Initial beta [deg]  = ', beta*180.0/PI

        write(*,'(a)') 'Newton Solver Settings:'
        write(*,'(a,1x,e25.16)') 'Finite Difference Step Size = ', finite_diff_step
        write(*,'(a,1x,e25.16)') '          Relaxation Factor = ', relax_factor
        write(*,'(a,1x,e25.16)') '                  Tolerance = ', newton_tol

        trim_array = trim_algorithm(this, this%init_state(9), newton_tol, p_wind, trim_type, is_trim_sideslip_angle, &
                trim_sideslip_angle, trim_bank_angle, trim_elevation_angle, trim_climb_angle, is_elevation_or_climb, &
                relax_factor, trim_azimuth_angle, finite_diff_step, is_trim_sct_load_factor, trim_load_factor)
        if (is_bank_or_beta_for_shss == "beta" .and. trim_type == "shss") then 
            this%init_state(10:13) = euler_to_quat([trim_array(2), trim_elevation_angle, trim_azimuth_angle])
            write(*,'(A12,1X,E22.14)') "theta[deg]", trim_elevation_angle*180.0/PI
            write(*,'(A12,1X,E22.14)') "phi[deg]", trim_array(2)*180.0/PI
            write(*,'(A12,1X,E22.14)') "psi[deg]", trim_azimuth_angle*180.0/PI
            write(*,'(A12,1X,E22.14)') "alpha[deg]", trim_array(1)*180.0/PI
            write(*,'(A12,1X,E22.14)') "beta[deg]", beta*180.0/PI
            write(*,'(A12,1X,E22.14)') "p[deg/s]", trim_array(3)*180.0/PI
            write(*,'(A12,1X,E22.14)') "q[deg/s]", trim_array(4)*180.0/PI
            write(*,'(A12,1X,E22.14)') "r[deg/s]", trim_array(5)*180.0/PI
            write(*,'(A12,1X,E22.14)') "da[deg]", trim_array(6)*180.0/PI
            write(*,'(A12,1X,E22.14)') "de[deg]", trim_array(7)*180.0/PI
            write(*,'(A12,1X,E22.14)') "dr[deg]", trim_array(8)*180.0/PI
            write(*,'(A12,1X,E22.14)') "tau", trim_array(9)
        else
            beta = trim_array(2)
            this%init_state(10:13) = euler_to_quat([trim_bank_angle, trim_elevation_angle, trim_azimuth_angle])
            write(*,'(A12,1X,E22.14)') "theta[deg]", trim_elevation_angle*180.0/PI
            write(*,'(A12,1X,E22.14)') "phi[deg]", trim_bank_angle*180.0/PI
            write(*,'(A12,1X,E22.14)') "psi[deg]", trim_azimuth_angle*180.0/PI
            write(*,'(A12,1X,E22.14)') "alpha[deg]", trim_array(1)*180.0/PI
            write(*,'(A12,1X,E22.14)') "beta[deg]", trim_array(2)*180.0/PI
            write(*,'(A12,1X,E22.14)') "p[deg/s]", trim_array(3)*180.0/PI
            write(*,'(A12,1X,E22.14)') "q[deg/s]", trim_array(4)*180.0/PI
            write(*,'(A12,1X,E22.14)') "r[deg/s]", trim_array(5)*180.0/PI
            write(*,'(A12,1X,E22.14)') "da[deg]", trim_array(6)*180.0/PI
            write(*,'(A12,1X,E22.14)') "de[deg]", trim_array(7)*180.0/PI
            write(*,'(A12,1X,E22.14)') "dr[deg]", trim_array(8)*180.0/PI
            write(*,'(A12,1X,E22.14)') "tau", trim_array(9)
        end if  
        alpha = trim_array(1)
        this%init_state(1) = this%init_V*cos(alpha)*cos(beta)
        this%init_state(2) = this%init_V*sin(beta)
        this%init_state(3) = this%init_V*sin(alpha)*cos(beta)
        this%init_state(4:6) = trim_array(3:5)
        this%controls(1:4) = trim_array(6:9) 

        load_factor_temp = 0.0
        ! get load factor from that 
        if (is_trim_sct_load_factor) then 
            call pseudo_aero(this, this%init_state)
            write(*,*) "Forces: ", this%FM(1), this%FM(2), this%FM(3)
            write(*,*) "Moments: ", this%FM(4), this%FM(5), this%FM(6)
            gravity_temp = gravity_English(-this%init_state(9))
            euler_temp = (/trim_bank_angle, trim_elevation_angle, trim_azimuth_angle/)
            qt_temp = euler_to_quat(euler_temp)
            xdot_temp = quat_dependent_to_base((/this%init_state(1),this%init_state(2),this%init_state(3)/)&
            , (/qt_temp(1), qt_temp(2), qt_temp(3), qt_temp(4)/))        
            v_t_temp = sqrt(xdot_temp(1)**2+xdot_temp(2)**2)
            a_c_temp = v_t_temp**2/(R_E_English - this%init_state(9))
            load_factor_temp = calc_load_factor(this, this%FM(1), this%FM(3), alpha, this%mass, gravity_temp, a_c_temp)
            write(*,*) "Load Factor: ", load_factor_temp
        end if 

    end subroutine init_to_trim

    subroutine mass_inertia(this) 
        implicit none 
        type(vehicle_t) :: this
        real :: gravity
        real :: det_I
        real :: I_tilde(3,3) 
        real :: a11, a22, a33, a12, a13, a21, a23, a31, a32
        real :: weight 
        this%I = 0.0
        this%Iinv = 0.0
        gravity = gravity_English(0.0)
        call jsonx_get(this%j_vehicle, "mass.weight[lbf]", weight)
        this%mass = weight/gravity
        call jsonx_get(this%j_vehicle, "mass.Ixx[slug-ft^2]", this%I(1,1))
        call jsonx_get(this%j_vehicle, "mass.Iyy[slug-ft^2]", this%I(2,2))
        call jsonx_get(this%j_vehicle, "mass.Izz[slug-ft^2]", this%I(3,3))
        if (this%type == "aircraft") then
            call jsonx_get(this%j_vehicle, "mass.Ixy[slug-ft^2]", this%I(1,2),0.0)
            call jsonx_get(this%j_vehicle, "mass.Ixz[slug-ft^2]", this%I(3,1),0.0)
            call jsonx_get(this%j_vehicle, "mass.Iyz[slug-ft^2]", this%I(3,2),0.0)
            call jsonx_get(this%j_vehicle, "mass.h[slug-ft^2/s]", this%h, 0.0, 3)
            this%hx = this%h(1)
            this%hy = this%h(2)
            this%hz = this%h(3)
            this%I(1,2) = -this%I(1,2)
            this%I(3,1) = -this%I(3,1)
            this%I(3,2) = -this%I(3,2)
            this%I(2,1) = this%I(1,2)
            this%I(1,3) = this%I(3,1)
            this%I(2,3) = this%I(3,2)

            this%Ixxb = this%I(1,1)
            this%Iyyb = this%I(2,2)
            this%Izzb = this%I(3,3)
            this%Ixyb = this%I(1,2)
            this%Ixzb = this%I(1,3)
            this%Iyzb = this%I(2,3)

            a11 = this%I(1,1)
            a22 = this%I(2,2)
            a33 = this%I(3,3)
            a12 = this%I(1,2)
            a13 = this%I(1,3)
            a21 = this%I(2,1)
            a23 = this%I(2,3)
            a31 = this%I(3,1)
            a32 = this%I(3,2)
            det_I = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 &
            - a13*a22*a31 - a12*a21*a33 - a11*a23*a32
            I_tilde(1,1) = a22*a33 - a23*a32
            I_tilde(1,2) = a13*a32 - a12*a33
            I_tilde(1,3) = a12*a23 - a13*a22
            I_tilde(2,1) = a23*a31 - a21*a33
            I_tilde(2,2) = a11*a33 - a13*a31
            I_tilde(2,3) = a13*a21 - a11*a23
            I_tilde(3,1) = a21*a32 - a22*a31
            I_tilde(3,2) = a12*a31 - a11*a32
            I_tilde(3,3) = a11*a22 - a12*a21
            this%Iinv = 1/(det_I)*I_tilde
            this%h_gyro = 0.0
            this%h_gyro(1,2) = -this%hz
            this%h_gyro(1,3) = this%hy
            this%h_gyro(2,1) = this%hz
            this%h_gyro(2,3) = -this%hx
            this%h_gyro(3,1) = -this%hy
            this%h_gyro(3,2) = this%hx
            this%hdot_gyro = 0.0
            ! write(*,*) "I", this%I 
            ! write(*,*) "Iinv", this%Iinv
            ! write(*,*) "Ixxb", this%Ixxb
            ! write(*,*) "Iyyb", this%Iyyb
            ! write(*,*) "Izzb", this%Izzb
            ! write(*,*) "Ixyb", this%Ixyb
            ! write(*,*) "Ixzb", this%Ixzb
            ! write(*,*) "Iyzb", this%Iyzb
        else if (this%type == "arrow") then
            this%I(1,2) = -this%I(1,2)
            this%I(3,1) = -this%I(3,1)
            this%I(3,2) = -this%I(3,2)
            this%I(2,1) = this%I(1,2)
            this%I(1,3) = this%I(3,1)
            this%I(2,3) = this%I(3,2)
            this%Ixxb = this%I(1,1)
            this%Iyyb = this%I(2,2)
            this%Izzb = this%I(3,3)
            this%Ixyb = this%I(1,2)
            this%Ixzb = this%I(1,3)
            this%Iyzb = this%I(2,3)
            this%Iinv(1,1) = 1/this%I(1,1)
            this%Iinv(2,2) = 1/this%I(2,2)
            this%Iinv(3,3) = 1/this%I(3,3)

            this%hx = 0.0
            this%hy = 0.0
            this%hz = 0.0
            this%h_gyro = 0.0
            this%hdot_gyro = 0.0 
        else if (this%type == "sphere") then
            if (abs(this%I(1,1))>0.000000001) then 
                this%Iinv(1,1) = 1.0/this%I(1,1)
            else 
                write(*,*) "this%I(1,1) is exactly zero, so the inverse is singular"
            end if 
            if (abs(this%I(2,2))>0.000000001) then 
                this%Iinv(2,2) = 1.0/this%I(2,2)
            else 
                write(*,*) "this%I(2,2) is exactly zero, so the inverse is singular"
            end if 
            if (abs(this%I(3,3))>0.000000001) then 
                this%Iinv(3,3) = 1.0/this%I(3,3)
            else 
                write(*,*) "this%I(3,3) is exactly zero, so the inverse is singular"
            end if 
            this%I(1,2) = -this%I(1,2)
            this%I(3,1) = -this%I(3,1)
            this%I(3,2) = -this%I(3,2)
            this%I(2,1) = this%I(1,2)
            this%I(1,3) = this%I(3,1)
            this%I(2,3) = this%I(3,2)
            this%Ixxb = this%I(1,1)
            this%Iyyb = this%I(2,2)
            this%Izzb = this%I(3,3)
            this%Ixyb = this%I(1,2)
            this%Ixzb = this%I(1,3)
            this%Iyzb = this%I(2,3)
            this%hx = 0.0
            this%hy = 0.0
            this%hz = 0.0
            this%h_gyro = 0.0
            this%hdot_gyro = 0.0 
        end if 
    end subroutine mass_inertia

    subroutine pseudo_aero(this, y) 
        implicit none
        type(vehicle_t) :: this
        real, intent(in) :: y(13)
        real :: da, de, dr, tau
        real :: V, alpha, beta, pbar, qbar, rbar, ahat
        real :: CL1, CL, CS, CD, Cll, Cm, Cn 
        real :: sa, ca, sb, cb, sign_a
        real :: Z, T, P, rho, a, mu
        real :: exp_pos, exp_neg, sigma, Reyn
        real :: CL_newt, CD_newt, Cm_newt
        real :: CL1_blended, CD1_blended, Cm1_blended
        real :: unit_V(3)
        real :: beta_f
        V = sqrt(y(1)**2 + y(2)**2 + y(3)**2)
        beta = asin(y(2)/V) ! Eq. 3.4.5
        alpha = atan2(y(3), y(1)) ! Eq. 3.4.4
        sa = sin(alpha)
        ca = cos(alpha)
        sb = sin(beta)
        cb = cos(beta)
        sign_a = sign(1.0, alpha)
        call std_atm_English(-y(9), Z, T, P, rho, a, mu)
        if (this%type == "aircraft") then 
            ahat = 0.0
            !!!! receive controls from python script here !!!!
            da = this%controls(1)
            de = this%controls(2)
            dr = this%controls(3)
            tau = this%controls(4)
            ! Clamp only negative throttle values
            ! if (tau < 0.0) then 
            !     tau = 0.0
            ! end if
            ! Allow tau > 1.0 for afterburner in trim calculations
            pbar = 0.5*y(4)*this%lat_ref/(V)
            qbar = 0.5*y(5)*this%long_ref/(V)
            rbar = 0.5*y(6)*this%lat_ref/(V)

            CL1 = this%CL0 +this%CLa*alpha
            CL = CL1 + this%CLqbar*qbar+this%CLahat*ahat + this%CLde*de
            CS = this%CSb*beta + (this%CSpbar+this%CSapbar*alpha)*pbar + this%CSrbar*rbar + this%CSda*da + this%CSdr*dr
            CD = this%CDL0 + this%CDL1*CL1 + this%CDL2*CL1**2 + this%CDS2*CS**2 + &
                (this%CDqbar + this%CDaqbar*alpha)*qbar + (this%CDde + this%CDade*alpha)*de + this%CDde2*de**2
            Cll = this%Clb*beta + this%Clpbar*pbar + (this%Clrbar + this%Clarbar*alpha)*rbar + this%Clda*da + this%Cldr*dr
            Cm = this%Cm0 + this%Cma*alpha + this%Cmqbar*qbar + this%Cmahat*ahat + this%Cmde*de 
            ! write(*,*) "Cm0", this%Cm0
            ! write(*,*) "Cma", this%Cma
            ! write(*,*) "alpha", alpha
            ! write(*,*) "Cmqbar", this%Cmqbar
            ! write(*,*) "qbar", qbar
            ! write(*,*) "Cmahat", this%Cmahat
            ! write(*,*) "ahat", ahat
            ! write(*,*) "Cmde", this%Cmde
            ! write(*,*) "de", de

            Cn = this%Cnb*beta + (this%Cnpbar + this%Cnapbar*alpha)*pbar + this%Cnrbar*rbar + &
                (this%Cnda + this%Cnada*alpha)*da + this%Cndr*dr
            if (this%include_stall) then
                ! CL 
                CL_newt = 2.0*sign_a*sa*sa*ca
                exp_pos = exp( this%lambda_CL*(alpha-this%alpha_0CL+this%alpha_sCL))
                exp_neg = exp(-this%lambda_CL*(alpha-this%alpha_0CL-this%alpha_sCL))
                sigma = (1.0 + exp_neg + exp_pos)/((1.0 + exp_neg)*(1.0 + exp_pos))
                CL = (1.0-sigma)*CL + sigma*CL_newt
                ! CD 
                CD_newt = 2.0*sin(abs(alpha))**3
                exp_pos = exp( this%lambda_CD*(alpha-this%alpha_0CD+this%alpha_sCD))
                exp_neg = exp(-this%lambda_CD*(alpha-this%alpha_0CD-this%alpha_sCD))
                sigma = (1.0 + exp_neg + exp_pos)/((1.0 + exp_neg)*(1.0 + exp_pos))
                CD = (1.0-sigma)*CD + sigma*CD_newt
                ! Cm
                Cm_newt = this%Cmmin*sign_a*sa*sa
                exp_pos = exp( this%lambda_Cm*(alpha-this%alpha_0Cm+this%alpha_sCm))
                exp_neg = exp(-this%lambda_Cm*(alpha-this%alpha_0Cm-this%alpha_sCm))
                sigma = (1.0 + exp_neg + exp_pos)/((1.0 + exp_neg)*(1.0 + exp_pos))
                Cm = (1.0-sigma)*Cm + sigma*Cm_newt
            end if 
            this%FM(1) = CL*sa-CS*ca*sb-CD*ca*cb
            this%FM(2) = CS*cb-CD*sb
            this%FM(3) = -CL*ca-CS*sa*sb-CD*sa*cb
            this%FM(4) = this%lat_ref*Cll
            this%FM(5) = this%long_ref*Cm
            this%FM(6) = this%lat_ref*Cn
            this%FM = 0.5*rho*V**2*this%sref*this%FM
            this%FM(1) = this%FM(1) + tau*this%Thrust0*(rho/this%rho0)**this%Ta
            this%FM(4:6) = this%FM(4:6) + cross_product_3D(this%aero_ref_location, this%FM(1:3)) ! body-fixed axis
        else if (this%type == "arrow") then 
            alpha = atan2(y(3), y(1)) ! Eq. 3.4.4
            beta = asin(y(2)/V) ! Eq. 3.4.5
            beta_f = atan2(y(2), y(1)) ! Eq. 3.4.13
            pbar = 0.5*y(4)*this%lat_ref/(V)
            qbar = 0.5*y(5)*this%long_ref/(V)
            rbar = 0.5*y(6)*this%lat_ref/(V)
            CL = this%CLa*alpha 
            CS = -this%CLa*beta_f 
            CD = this%CDL0 + this%CDL2*CL**2 + this%CDL2*CS**2
            Cll = this%Cll0 + this%Clpbar*pbar
            Cm = this%Cma*alpha + this%Cmqbar*qbar 
            Cn = -this%Cma*beta_f + this%Cmqbar*rbar 
            ! Eq. 5.2.3 in book. 
            this%FM(1) = CL*sin(alpha)-CS*cos(alpha)*sin(beta)-CD*cos(alpha)*cos(beta)
            this%FM(2) = CS*cos(beta)-CD*sin(beta)
            this%FM(3) = -CL*cos(alpha)-CS*sin(alpha)*sin(beta)-CD*sin(alpha)*cos(beta)
            ! Eq. 5.2.4 in book
            this%FM(4) = this%lat_ref*Cll
            this%FM(5) = this%long_ref*Cm
            this%FM(6) = this%lat_ref*Cn
            this%FM = 0.5*rho*V**2*this%sref*this%FM
            this%FM(1) = this%FM(1) + tau*this%Thrust0*(rho/this%rho0)**this%Ta
            this%FM(4:6) = this%FM(4:6) + cross_product_3D(this%aero_ref_location, this%FM(1:3)) ! body-fixed axis
        else if (this%type == "sphere") then
            this%FM(4:6) = 0.0
            Reyn = 2*rho*V*this%long_ref/mu
            unit_V(1:3) = y(1:3)
            unit_V = unit_V/V
            if (Reyn < 0.01) then 
                CD = 2405.0
            else if (0.01 <= Reyn .and. Reyn <= 450000.0) then 
                CD = 24/Reyn + 6/(1+sqrt(Reyn)) + 0.4
            else if (450000.0 <= Reyn .and. Reyn <= 560000.0) then 
                CD = 1.0*10.0**29*Reyn**(-5.211)
            else if (560000.0 <= Reyn .and. Reyn <= 14000000.0) then 
                CD = -2.0*10.0**(-23)*Reyn**3 - 1.0*10.0**(-16)*Reyn**2 + 9.0*10.0**(-9)*Reyn + 0.069
            else 
                CD = 0.12
            end if 
            this%FM(1:3) = -0.5*rho*V**2*PI*this%long_ref**2*CD*unit_V
        else
            write(*,*) "NO VALID VEHICLE TYPE FOUND IN JSON FILE. THIS ERROR MESSAGE IS IN THE pseudo_aero subroutine!"
        end if 
        if (this%rk4_verbose) then
            write(*,*) "| pseudo aerodynamics (F,M) = "
            write(*,'(14E22.13)') this%FM
        end if
    end subroutine pseudo_aero

    function differential_equations(this, time, state) result(res)
        implicit none 
        type(vehicle_t) :: this
        real, intent(in) :: time
        real, intent(in), dimension(:) :: state
        real :: res(size(state))
        real :: u,v,w,p,q,r,x,y,z
        real :: e0,ex,ey,ez
        real, dimension(3) :: pqr_temp, rot_and_inertia_temp
        real :: gravity, a_c, v_t
        u = state(1)
        v = state(2)
        w = state(3)
        p = state(4)
        q = state(5)
        r = state(6)
        pqr_temp = [p,q,r]
        x = state(7)
        y = state(8)
        z = state(9)
        e0 = state(10)
        ex = state(11)
        ey = state(12)
        ez = state(13)
        gravity = gravity_English(-z)
        call pseudo_aero(this, state)
        rot_and_inertia_temp = & 
        [this%FM(4) + dot_product(this%h_gyro(1,:),pqr_temp) + ((this%Iyyb-this%Izzb)*q*r-this%Iyzb*(q**2-r**2)-&
        this%Ixzb*p*q+this%Ixyb*p*r)-this%hdot_gyro(1),&
        this%FM(5) + dot_product(this%h_gyro(2,:),pqr_temp) + ((this%Izzb-this%Ixxb)*p*r-this%Ixzb*(r**2-p**2)-&
        this%Ixyb*q*r+this%Iyzb*p*q)-this%hdot_gyro(2),& 
        this%FM(6) + dot_product(this%h_gyro(3,:),pqr_temp) + ((this%Ixxb-this%Iyyb)*p*q-this%Ixyb*(p**2-q**2)-&
        this%Iyzb*p*r+this%Ixzb*q*r)-this%hdot_gyro(3)]
        ! get xdot ydot zdot first so we can add gravity relief to udot, vdot, wdot
        res(7:9) = quat_dependent_to_base((/u,v,w/), (/e0, ex, ey, ez/)) !!! + wind ! xdot ydot zdot earth-fixed
        !!!! velocity with gravity relief
        v_t = sqrt(res(7)**2+res(8)**2)
        a_c = v_t**2/(R_E_English - z)
        !!!!!
        res(1) = 1/this%mass * this%FM(1) + (gravity - a_c) * 2 *(ex*ez-ey*e0) + r*v - q*w ! udot body-fixed
        res(2) = 1/this%mass * this%FM(2) + (gravity - a_c) * 2 *(ey*ez+ex*e0) + p*w - r*u ! vdot body-fixed
        res(3) = 1/this%mass * this%FM(3) + (gravity - a_c) * (ez**2+e0**2-ex**2-ey**2) + q*u - p*v ! wdot body-fixed
        res(4:6) = matmul(this%Iinv, rot_and_inertia_temp) ! pdot, qdot, rdot body-fixed
        res(10) = 0.5 * dot_product((/-ex, -ey, -ez/),pqr_temp) !e0
        res(11) = 0.5 * dot_product((/e0, -ez, ey/),pqr_temp) !ex
        res(12) = 0.5 * dot_product((/ez, e0, -ex/),pqr_temp) !ey
        res(13) = 0.5 * dot_product((/-ey, ex, e0/),pqr_temp) !ez
        if (this%rk4_verbose) then
            write(*,*) " | diff eq results = "
            write(*,'(14E22.13)') res
            ! write(*,*) ""
            ! write(*,*) "pqr", pqr_temp
            ! write(*,*) "h_gyro", this%h_gyro
            ! write(*,*) "hdot_gyro", this%hdot_gyro
            ! write(*,*) "I", this%I
            ! write(*,*) "Iinv", this%Iinv
            ! write(*,*) "Ixxb", this%Ixxb
            ! write(*,*) "Iyyb", this%Iyyb
            ! write(*,*) "Izzb", this%Izzb
            ! write(*,*) "Iyzb", this%Iyzb
            ! write(*,*) "Ixzb", this%Ixzb
            ! write(*,*) "Ixyb", this%Ixyb
            ! write(*,*) "rot_and_inertia_temp", rot_and_inertia_temp
            write(*,*) ""
        end if
    end function differential_equations

    function cross_product_3D(a, b) result(result)
        implicit none
        real, intent(in) :: a(3), b(3) 
        real :: result(3)
        result(1) = a(2)*b(3) - a(3)*b(2)
        result(2) = a(3)*b(1) - a(1)*b(3)
        result(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product_3D

    subroutine vehicle_tick_state(this, t,dt)
        implicit none 
        type(vehicle_t) :: this
        real, intent(in) ::  t, dt
        real :: y(13), y1(13)
        real :: euler_angles(3), azimuth 
        real :: print_times(3)

        print_times = (/6.63,6.64,10.0/)
        y = this%state 
        y1 = runge_kutta(this,t,y,dt)
        if (sqrt(y1(4)**2+y1(5)**2+y1(6)**2)/2.0/PI*dt>0.1) write(*,*) &
        'Warning, rotation rates large for Rk4. See Eq. 5.7.3 in the book'
        if (geographic_model_ID>0) call update_geographic(this,y,y1)
        call quat_norm(y1(10:13))
        this%state = y1 
        if (this%save_states .and. geographic_model_ID>0 .and. &
        (abs(t+dt-print_times(1))<0.001 .or. (abs(t+dt-print_times(2))<0.001) &
        .or. (abs(t+dt-print_times(3))<0.001))) then
        ! if (this%save_states .and. geographic_model_ID>0) then
            euler_angles = quat_to_euler(this%state(10:13))
            azimuth = euler_angles(3)
            write(this%iunit_geographic,'(*(G0.8,:,","))') t+dt, this%latitude*180.0/PI, this%longitude*180.0/PI, azimuth*180.0/PI
        end if
         if (this%save_states .and. &
        (abs(t+dt-print_times(1))<0.001 .or. (abs(t+dt-print_times(2))<0.001) &
        .or. (abs(t+dt-print_times(3))<0.001))) then
        !  if (this%save_states) then 
            call vehicle_write_state(this, t+dt,y1)
        end if
    end subroutine vehicle_tick_state

    subroutine update_geographic(this, y1, y2)
        implicit none 
        type(vehicle_t) :: this 
        real, intent(in) :: y1(13)
        real, intent(inout) :: y2(13)
        real :: dx, dy, dz, d 
        real :: theta, g1, xhat, yhat, zhat, xhp, yhp, zhp, rhat, Chat, Shat 
        real :: Phi1, Psi1, H1 
        real :: cP, cT, sP, sT, cs, cg, sg, dg
        real :: quat(4) 
        real :: Rp, Re, e2 
        real :: temp, Rx, Ry, tx, ty

        ! write(*,*) "state 1"
        ! write(*,*) y1 
        ! write(*,*) "state 2 in"
        ! write(*,*) y2 

        dx = y2(7) - y1(7)
        dy = y2(8) - y1(8)
        dz = y2(9) - y1(9)

        d = sqrt(dx**2 + dy**2)
        if(d<1e-14) then 
            ! don't do anything
        else 
            H1 = -y1(9)
            Phi1 = this%latitude 
            Psi1 = this%longitude 
            cP = cos(Phi1)
            sP = sin(Phi1)
            if (geographic_model_ID == 1) then ! spherical model
                theta = d/(R_E/0.3048 + H1 - 0.5*dz)
                cT = cos(theta)
                sT = sin(theta)
                g1 = atan2(dy,dx)
                cg = cos(g1) 
                sg = sin(g1)
                xhat = cP*cT - sP*sT*cg 
                yhat = sT*sg 
                zhat = sP*cT + cP*sT*cg
                xhp = -cP*sT - sP*cT*cg 
                yhp = cT*sg 
                zhp = -sP*sT +cP*cT*cg 
                rhat = sqrt(xhat**2 + yhat**2)
                ! write(*,*) "Latitude before", this%latitude 
                ! write(*,*) "Longitude before", this%longitude 
                this%latitude = atan2(zhat,rhat)
                this%longitude = Psi1 + atan2(yhat, xhat)
                Chat = xhat**2*zhp 
                Shat = (xhat*yhp - yhat*xhp)*cos(this%latitude)**2*cos(this%longitude-Psi1)**2
                dg = atan2(Shat,Chat) - g1
                ! write(*,*) "H1", H1
                ! write(*,*) "Phi1", Phi1
                ! write(*,*) "Psi1", Psi1
                ! write(*,*) "cP", cP
                ! write(*,*) "sP", sP
                ! write(*,*) "cT", cT
                ! write(*,*) "sT", sT
                ! write(*,*) "g1", g1
                ! write(*,*) "cg", cg
                ! write(*,*) "sg", sg
                ! write(*,*) "xhat", xhat 
                ! write(*,*) "yhat", yhat 
                ! write(*,*) "zhat", zhat
                ! write(*,*) "xhp", xhp  
                ! write(*,*) "yhp", yhp  
                ! write(*,*) "zhp", zhp  
                ! write(*,*) "rhat", rhat
                ! write(*,*) "Chat", Chat 
                ! write(*,*) "Shat", Shat
                ! write(*,*) "dg", dg


            else ! ellipse  
                Rp = 6356.7516/0.3048*1000.0
                Re = 6378.1363/0.3048*1000.0
                e2 = 1.0 - (Rp/Re)**2 
                temp = 1.0 - e2*sin(Phi1)**2 
                Rx = Re*(1.0-e2)/(temp**1.5)
                Ry = Re/(temp**0.5)
                tx = dx/(Rx + H1 - 0.5*dz)
                ty = dy/(Ry + H1 - 0.5*dz)
                xhat = (1.0-e2)*(cos(Phi1 + tx) - cos(Phi1)) + temp*cos(ty)*cos(Phi1)
                yhat = temp*sin(ty)
                zhat = (1.0-e2)*(sin(Phi1 + tx) - sin(Phi1)) + temp*(cos(ty)-e2)*sin(Phi1)
                rhat = sqrt(xhat**2 + yhat**2) 
                this%latitude = atan2(zhat, (1.0-e2)*rhat)
                this%longitude = Psi1 + atan2(yhat, xhat) 

                dg = (this%longitude - Psi1)*sin(0.5*(this%latitude + Phi1))*(1.0 - e2)/temp
            end if 
            
            if(this%longitude > PI) this%longitude = this%longitude - 2.0*PI 
            if(this%longitude < -PI) this%longitude = this%longitude + 2.0*PI 

            ! write(*,*) "Latitude after", this%latitude 
            ! write(*,*) "Longitude after", this%longitude 

            cg = cos(0.5*dg)
            sg = sin(0.5*dg)
            
            ! write(*,*) "cg second", cg
            ! write(*,*) "sg second", sg

            quat(1) = -y2(13)
            quat(2) = -y2(12)
            quat(3) =  y2(11)
            quat(4) =  y2(10)
            y2(10:13) = cg*y2(10:13) + sg*quat(:)
            
            ! write(*,*)
            ! write(*,*) "state 2 out"
            ! write(*,*) y2 

        end if 
    end subroutine update_geographic

    subroutine vehicle_write_state(this, time, state)
        implicit none 
        type(vehicle_t) :: this
        real, intent(in) :: time, state(13)
        logical :: is_open
        ! inquire(file=this%states_filename, opened = is_open)
        ! if (is_open) then 
            ! write(*,*) this%states_filename, ' is already open. THIS PROGRAM WILL NOT RUN IF THAT IS NOT CLOSED'
        ! else 
        write(this%iunit_states,'(*(G0.15,:,","))') time, state
        ! end if 
    end subroutine vehicle_write_state

    function runge_kutta(this, t_0, state, delta_t) result(state_out)
        implicit none 
        type(vehicle_t) :: this
        real, intent(in) :: t_0
        real, intent(in), dimension(:) :: state
        real, intent(in) :: delta_t
        real :: state_out(size(state))
        real :: k1(13), k2(13), k3(13), k4(13)
        real :: state_temp(13)
        real :: t_0_plus_delta_t_over_2, delta_t_over_2
        delta_t_over_2 = delta_t/2.0
        if (this%rk4_verbose) then
            write(*,*) "----------------BEGINNING OF SINGLE RK4 INTEGRATION CALL---------------------"
            write(*,*) "state of the vehicle at the beginning of this RK4 integration step:"
            write(*,*) "t   u   v   w   p   q   r   x   y   z   e0   ex   ey   ez"
            write(*,'(14E22.13)') t_0, state 
            write(*,*) ""
            write(*,*) "RK4 Function called..."
            write(*,*) ""
            t_0_plus_delta_t_over_2 = t_0 + delta_t_over_2
            write(*,*) "diff_eq function called... " 
            write(*,*) "  | RK4 call number        =         1"
            write(*,*) "  | time [s]               =          ", t_0
            write(*,*) "  | State Vector Coming in =          "
            write(*,'(14E22.13)') state
            k1 = differential_equations(this,t_0, state)
                state_temp = state
                state_temp = state_temp + delta_t_over_2 * k1
            write(*,*) "diff_eq function called... " 
            write(*,*) "  | RK4 call number        =         2"
            write(*,*) "  | time [s]               =          ", t_0_plus_delta_t_over_2
            write(*,*) "  | State Vector Coming in =          "
            write(*,'(14E22.13)') state_temp
            k2 = differential_equations(this,t_0_plus_delta_t_over_2, state_temp)
                state_temp = state
                state_temp = state_temp + delta_t_over_2 * k2
            write(*,*) "diff_eq function called... " 
            write(*,*) "  | RK4 call number        =         3"
            write(*,*) "  | time [s]               =          ", t_0_plus_delta_t_over_2
            write(*,*) "  | State Vector Coming in =          "
            write(*,'(14E22.13)') state_temp
            k3 = differential_equations(this,t_0_plus_delta_t_over_2, state_temp)
                state_temp = state 
                state_temp = state_temp + delta_t * k3 
            write(*,*) "diff_eq function called... " 
            write(*,*) "  | RK4 call number        =         4"
            write(*,*) "  | time [s]               =          ", t_0 + delta_t 
            write(*,*) "  | State Vector Coming in =          "
            write(*,'(14E22.13)') state_temp
            k4 = differential_equations(this,t_0 + delta_t, state_temp)
            state_out = state + (delta_t/6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
            write(*,*) "state of the vehicle after this RK4 integration step:"
            write(*,*) "t   u   v   w   p   q   r   x   y   z   e0   ex   ey   ez"
            write(*,'(14E22.13)') state_out
            write(*,*) "----------------END OF SINGLE RK4 INTEGRATION CALL---------------------"
        else
            t_0_plus_delta_t_over_2 = t_0 + delta_t_over_2
            k1 = differential_equations(this,t_0, state)
                state_temp = state
                state_temp = state_temp + delta_t_over_2 * k1
            k2 = differential_equations(this,t_0_plus_delta_t_over_2, state_temp)
                state_temp = state
                state_temp = state_temp + delta_t_over_2 * k2
            k3 = differential_equations(this,t_0_plus_delta_t_over_2, state_temp)
                state_temp = state 
                state_temp = state_temp + delta_t * k3 
            k4 = differential_equations(this,t_0 + delta_t, state_temp)
            state_out = state + (delta_t/6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        end if
    end function runge_kutta

    !!!!!! trim stuff !!!!!!
    function trim_algorithm(this, H_altitude, newton_tol, p_wind, trim_type, is_trim_sideslip_angle, &
        trim_sideslip_angle, trim_bank_angle, trim_elevation_angle, trim_climb_angle, &
        is_elevation_or_climb, relax_factor, trim_azimuth_angle, finite_diff_step, is_trim_sct_load_factor, &
        trim_load_factor) result(trim_result)
        implicit none 
        type(vehicle_t) :: this
        real, intent(in) :: H_altitude, newton_tol
        real, intent(in) :: p_wind

        character(len=*), intent(in) :: trim_type, is_elevation_or_climb
        logical, intent(in) :: is_trim_sideslip_angle, is_trim_sct_load_factor
        real, intent(in) :: trim_sideslip_angle, trim_climb_angle
        real, intent(in) :: relax_factor, trim_azimuth_angle, finite_diff_step
        real, intent(in) :: trim_load_factor
        real :: trim_elevation_angle, trim_bank_angle
        real :: trim_result(9)
        real :: alpha, beta, p, q, r, da, de, dr, tau,x,y,z
        real :: pos(3), quat_orientation(4)
        real :: phi, theta, psi
        real :: current_error
        real :: u,v,w, sct_pqr_coeff
        real, allocatable :: DeltaG(:)
        real :: residual(6)
        real :: jacobian(6,6)
        real :: gravity, a_c, v_t 
        real :: newton_input(6)
        integer :: i, j
        real :: load_factor_error, p_old, q_old, r_old
        real :: phi_lf_convergence_tol
        integer :: k_lf
        real :: xdot_temp(3), euler_temp(3), qt(4)
        write(*,*) ""
        write(*,*) "Beginning trim algorithm..."
        gravity = gravity_English(-H_altitude)
        alpha = 0.0
        allocate(DeltaG(6))
        if (trim_type == "shss" .and. is_trim_sideslip_angle) then 
            beta = trim_sideslip_angle
        else
            beta = 0.0
        end if 
        u = this%init_V*cos(alpha)*cos(beta)
        v = this%init_V*sin(beta)
        w = this%init_V*sin(alpha)*cos(beta) 
        x = 0.0
        y = 0.0
        z = H_altitude
        p = 0.0
        q = 0.0
        r = 0.0
        da = 0.0
        de = 0.0 
        dr = 0.0
        tau = 0.0
        if (trim_type == "shss") then
            if (is_trim_sideslip_angle) then 
                beta = trim_sideslip_angle
                phi = trim_bank_angle
            else 
                phi = trim_bank_angle
            end if 
        else
            phi = trim_bank_angle
        end if 
        if (trim_type == "sct") then 
            if (is_elevation_or_climb == "climb") then 
                trim_elevation_angle = calc_theta_from_climb_angle(this, u, v, w, phi, trim_climb_angle)
            end if
            theta = trim_elevation_angle
            psi = trim_azimuth_angle
            ! Converge phi and pqr when load factor trim is active
            if (is_trim_sct_load_factor) then
                phi_lf_convergence_tol = 1.0e-9
                load_factor_error = 100.0
                k_lf = 1
                do while (load_factor_error > phi_lf_convergence_tol .and. k_lf < 100)
                    ! Calculate p, q, r based on current phi
                    euler_temp = (/phi, theta, psi/)
                    qt = euler_to_quat(euler_temp)
                    xdot_temp = quat_dependent_to_base((/u,v,w/), (/qt(1), qt(2), qt(3), qt(4)/))        
                    v_t = sqrt(xdot_temp(1)**2+xdot_temp(2)**2)
                    a_c = v_t**2/(R_E_English - z)
                    sct_pqr_coeff = (gravity-a_c)*sin(phi)*cos(theta)/(u*cos(theta)*cos(phi)+w*sin(theta))
                    p_old = p
                    q_old = q
                    r_old = r
                    p = -sct_pqr_coeff*(sin(theta))
                    q = sct_pqr_coeff*(sin(phi)*cos(theta))
                    r = sct_pqr_coeff*(cos(phi)*cos(theta))
                    ! Recalculate phi from updated p, q, r
                    trim_bank_angle = calc_phi_from_load_factor(this, u,v,w,p,q,r,alpha,theta, gravity, trim_load_factor, z)
                    phi = trim_bank_angle
                    ! Check convergence
                    load_factor_error = max(abs(p-p_old), abs(q-q_old), abs(r-r_old))
                    k_lf = k_lf + 1
                end do
            else
                ! No load factor constraint, just calculate pqr once
                euler_temp = (/phi, theta, psi/)
                qt = euler_to_quat(euler_temp)
                xdot_temp = quat_dependent_to_base((/u,v,w/), (/qt(1), qt(2), qt(3), qt(4)/))        
                v_t = sqrt(xdot_temp(1)**2+xdot_temp(2)**2)
                a_c = v_t**2/(R_E_English - z)
                ! Check for potential divide-by-zero
                if (abs(u*cos(theta)*cos(phi)+w*sin(theta)) < 1.0e-10) then
                    write(*,*) "WARNING: Near-zero denominator in sct_pqr_coeff calculation!"
                    write(*,*) "  u=", u, " w=", w, " theta[deg]=", theta*180.0/PI, " phi[deg]=", phi*180.0/PI
                    write(*,*) "  Denominator=", u*cos(theta)*cos(phi)+w*sin(theta)
                end if
                ! Check for potential divide-by-zero
                if (abs(u*cos(theta)*cos(phi)+w*sin(theta)) < 1.0e-10) then
                    write(*,*) "WARNING: Near-zero denominator in sct_pqr_coeff calculation!"
                    write(*,*) "  u=", u, " w=", w, " theta[deg]=", theta*180.0/PI, " phi[deg]=", phi*180.0/PI
                    write(*,*) "  Denominator=", u*cos(theta)*cos(phi)+w*sin(theta)
                end if
                sct_pqr_coeff = (gravity-a_c)*sin(phi)*cos(theta)/(u*cos(theta)*cos(phi)+w*sin(theta))
                p = -sct_pqr_coeff*(sin(theta))
                q = sct_pqr_coeff*(sin(phi)*cos(theta))
                r = sct_pqr_coeff*(cos(phi)*cos(theta))
            end if
        else if (trim_type == "vbr") then 
            p = (p_wind/this%init_V)*u
            q = (p_wind/this%init_V)*v
            r = (p_wind/this%init_V)*w
        else
            p = 0.0
            q = 0.0
            r = 0.0
        end if 
        current_error = 100.0
        if (trim_type == "shss" .and. is_trim_sideslip_angle) then
                newton_input = [alpha, phi, da, de, dr, tau]
                write(*,*) "newton initial input", newton_input
            else
                newton_input = [alpha, beta, da, de, dr, tau]
                write(*,*) "newton initial input", newton_input
        end if
        j = 1
        do while(current_error > newton_tol)
            ! Temporarily disable throttle clamping to diagnose convergence issue
            ! if (newton_input(6) < 0.0) then 
            !     newton_input(6) = 0.0
            ! else if (newton_input(6) > 1.0) then 
            !     newton_input(6) = 1.0
            ! end if 
            alpha = newton_input(1)
            if (trim_type == "shss" .and. is_trim_sideslip_angle) then 
                beta = trim_sideslip_angle
                phi = newton_input(2)
            else
                beta = newton_input(2)
            end if 
            ! beta = newton_input(2)
            u = this%init_V*cos(alpha)*cos(beta)
            v = this%init_V*sin(beta)
            w = this%init_V*sin(alpha)*cos(beta) 
            if (is_elevation_or_climb == "climb") then 
                trim_elevation_angle = calc_theta_from_climb_angle(this, u, v, w, phi, trim_climb_angle)
            end if 
            ! write(*,*) "Trim theta", trim_elevation_angle * 180.0/PI this is in degrees
            if (trim_type == "sct") then  
                theta = trim_elevation_angle
                psi = trim_azimuth_angle
                ! Converge phi and pqr when load factor trim is active
                if (is_trim_sct_load_factor) then
                    phi_lf_convergence_tol = 1.0e-9
                    load_factor_error = 100.0
                    k_lf = 1
                    do while (load_factor_error > phi_lf_convergence_tol .and. k_lf < 100)
                        ! Calculate p, q, r based on current phi
                        euler_temp = (/phi, theta, psi/)
                        qt = euler_to_quat(euler_temp)
                        xdot_temp = quat_dependent_to_base((/u,v,w/), (/qt(1), qt(2), qt(3), qt(4)/))        
                        v_t = sqrt(xdot_temp(1)**2+xdot_temp(2)**2)
                        a_c = v_t**2/(R_E_English - z)
                        sct_pqr_coeff = (gravity-a_c)*sin(phi)*cos(theta)/(u*cos(theta)*cos(phi)+w*sin(theta))
                        p_old = p
                        q_old = q
                        r_old = r
                        p = -sct_pqr_coeff*(sin(theta))
                        q = sct_pqr_coeff*(sin(phi)*cos(theta))
                        r = sct_pqr_coeff*(cos(phi)*cos(theta))
                        ! Recalculate phi from updated p, q, r
                        trim_bank_angle = calc_phi_from_load_factor(this, u,v,w,p,q,r,alpha,theta, gravity, trim_load_factor, z)
                        phi = trim_bank_angle
                        ! Check convergence
                        load_factor_error = max(abs(p-p_old), abs(q-q_old), abs(r-r_old))
                        k_lf = k_lf + 1
                    end do
                else
                    ! No load factor constraint, just calculate pqr once
                    euler_temp = (/phi, theta, psi/)
                    qt = euler_to_quat(euler_temp)
                    xdot_temp = quat_dependent_to_base((/u,v,w/), (/qt(1), qt(2), qt(3), qt(4)/))        
                    v_t = sqrt(xdot_temp(1)**2+xdot_temp(2)**2)
                    a_c = v_t**2/(R_E_English - z)
                    sct_pqr_coeff = (gravity-a_c)*sin(phi)*cos(theta)/(u*cos(theta)*cos(phi)+w*sin(theta))
                    p = -sct_pqr_coeff*(sin(theta))
                    q = sct_pqr_coeff*(sin(phi)*cos(theta))
                    r = sct_pqr_coeff*(cos(phi)*cos(theta))
                end if
            else if (trim_type == "vbr") then 
                p = (p_wind/this%init_V)*u
                q = (p_wind/this%init_V)*v
                r = (p_wind/this%init_V)*w
            else 
                p = 0.0
                q = 0.0
                r = 0.0
            end if 
            residual = calc_residual(this, newton_input, p, q, r, is_trim_sideslip_angle, &
                    trim_azimuth_angle, trim_bank_angle, trim_elevation_angle, trim_sideslip_angle, trim_type)
            jacobian = create_jacobian(this, newton_input, finite_diff_step, p, q, r, is_trim_sideslip_angle, &
                    trim_azimuth_angle, trim_bank_angle, trim_elevation_angle, trim_sideslip_angle, trim_type)
            call lu_solve(6,jacobian,residual,DeltaG)
            newton_input = newton_input - relax_factor*DeltaG
            residual = calc_residual(this, newton_input, p, q, r, is_trim_sideslip_angle, &
                        trim_azimuth_angle, trim_bank_angle, trim_elevation_angle, trim_sideslip_angle, trim_type)
            current_error = maxval(abs(residual))
            ! give diagnostic every 10 iterations and then every 100 iterations after 100 iterations
            if ((j <= 100 .and. mod(j,10) == 0) .or. (j > 100 .and. mod(j,100) == 0)) then
                write(*,'(A,I5,A,E12.5,A,6E12.4)') "Iter ", j, " Error: ", current_error, &
                    " State: ", newton_input(1)*180./PI, newton_input(2)*180./PI, newton_input(3:6)
            end if
            j = j + 1
            if (j > 1000) then
                write(*,*) "WARNING: Exceeded 1000 iterations without convergence!"
                write(*,*) "Current error:", current_error, " Tolerance:", newton_tol
                exit
            end if
        end do 
        write(*,*) "Trim converged!"
        write(*,*) "Converged in ", j-1, " iterations with error ", current_error
        alpha = newton_input(1) 
        beta = newton_input(2)
        da = newton_input(3)
        de = newton_input(4)
        dr = newton_input(5)
        tau = newton_input(6)
        u = this%init_V*cos(alpha)*cos(beta)
        v = this%init_V*sin(beta)
        w = this%init_V*sin(alpha)*cos(beta) 
        if (trim_type == "shss" .and. is_trim_sideslip_angle) then
            phi = newton_input(2)
        else 
            phi = trim_bank_angle
        end if
        theta = trim_elevation_angle
        psi = trim_azimuth_angle
        write(*,*) "Final Trim Results:"
        write(*,*) "Trim type ", trim_type
        if (trim_type == "sct") then
            euler_temp = (/phi, theta, psi/)
            qt = euler_to_quat(euler_temp)
            xdot_temp = quat_dependent_to_base((/u,v,w/), (/qt(1), qt(2), qt(3), qt(4)/))
            !!!! velocity with gravity relief
            v_t = sqrt(xdot_temp(1)**2+xdot_temp(2)**2)
            a_c = v_t**2/(R_E_English - z)
            sct_pqr_coeff = (gravity-a_c)*sin(phi)*cos(theta)/(u*cos(theta)*cos(phi)+w*sin(theta)) !!!!! INCLUDE GRAVITY RELIEF HERE
            p = -sct_pqr_coeff*(sin(theta))
            q = sct_pqr_coeff*(sin(phi)*cos(theta))
            r = sct_pqr_coeff*(cos(phi)*cos(theta))
        else if (trim_type == "vbr") then 
            p = (p_wind/this%init_V)*u
            q = (p_wind/this%init_V)*v
            r = (p_wind/this%init_V)*w
        else 
            p = 0.0
            q = 0.0
            r = 0.0
        end if 
        trim_result = [alpha, beta, p, q, r, da, de, dr, tau]
    end function trim_algorithm

    function calc_theta_from_climb_angle(this, u, v, w, phi, trim_climb_angle) result(return_theta)
        implicit none 
        type(vehicle_t) :: this
        real, intent(in) :: u, v, w, phi, trim_climb_angle
        real :: return_theta, error_for_elevation, temp_lhs, temp_rhs_1, temp_rhs_2
        real :: parenth_temp, S_theta_1, S_theta_2, theta_1, theta_2
        ! write(*,*) "Solving for elevation angle given a climb angle"
        error_for_elevation = 1e-12
        temp_lhs = this%init_V*sin(trim_climb_angle)
        parenth_temp = v*sin(phi)+w*cos(phi)
        ! Check for negative value under square root
        if (u**2 + parenth_temp**2 - this%init_V**2*sin(trim_climb_angle)**2 < 0.0) then
            write(*,*) "ERROR: Negative value under square root in calc_theta_from_climb_angle!"
            write(*,*) "  u=", u, " parenth_temp=", parenth_temp
            write(*,*) "  init_V=", this%init_V, " climb_angle[deg]=", trim_climb_angle*180.0/PI
            write(*,*) "  Value under sqrt=", u**2 + parenth_temp**2 - this%init_V**2*sin(trim_climb_angle)**2
        end if
        S_theta_1 = (u*this%init_V*sin(trim_climb_angle)+parenth_temp*sqrt(u**2 + parenth_temp**2 - &
        this%init_V**2*sin(trim_climb_angle)**2))/(u**2 + parenth_temp**2)
        theta_1 = asin(S_theta_1)
        ! write(*,*) "theta 1 [deg] = ", theta_1 * 180.0/PI
        S_theta_2 = (u*this%init_V*sin(trim_climb_angle)-parenth_temp*sqrt(u**2 + parenth_temp**2 - &
        this%init_V**2*sin(trim_climb_angle)**2))/(u**2 + parenth_temp**2)
        theta_2 = asin(S_theta_2)
        ! write(*,*) "theta 2 [deg] = ", theta_2 * 180.0/PI
        temp_rhs_1 = u*S_theta_1 - parenth_temp*cos(theta_1)
        temp_rhs_2 = u*S_theta_2 - parenth_temp*cos(theta_2)
        if (abs(temp_lhs-temp_rhs_1) <= error_for_elevation) then
            return_theta = theta_1
            ! write(*,*) "correct theta [deg] = ", theta_1 * 180.0/PI
        else if (abs(temp_lhs - temp_rhs_2) <= error_for_elevation) then 
            return_theta = theta_2 
                ! write(*,*) "correct theta [deg] = ", theta_2 * 180.0/PI
        else 
            write(*,*) "Error", abs(temp_lhs - temp_rhs_2)
            write(*,*) "WARNING, BOTH THETA VALUES DO NOT SATISFY THE LHS OF EQ. 7.2.9" 
        end if 
    end function calc_theta_from_climb_angle

    function calc_phi_from_load_factor(this, u, v, w, p, q, r, alpha, theta, g, load_factor, z) result(return_phi)
        implicit none 
        type(vehicle_t) :: this
        real, intent(in) :: u, v, w, p, q, r, alpha, theta, g, load_factor, z
        real :: C_theta, S_theta, C_alpha, S_alpha, A, error_tol, current_error
        real :: phi_guess, phi_current, g_minus_ac 
        real :: psi
        real :: return_phi, parenthesis_term
        integer :: i 
        real :: xdot_temp(3), euler_temp(3), qt(4)
        real :: a_c, v_t 
        current_error = 100.0
        error_tol = 1e-13
        a_c = 0.0
        g_minus_ac = g-a_c 
        C_theta = cos(theta)
        S_theta = sin(theta) 
        C_alpha = cos(alpha)
        S_alpha = sin(alpha) 
        parenthesis_term = (C_theta)/(load_factor)
        phi_guess = acos( max(-1.0, min(1.0, parenthesis_term)))
        i = 1
        psi = 0.0
        do while(current_error > error_tol .and. i < 500)
            euler_temp = (/phi_guess, theta, psi/)
            qt = euler_to_quat(euler_temp)
            xdot_temp = quat_dependent_to_base((/u,v,w/), (/qt(1), qt(2), qt(3), qt(4)/))
            !!!! velocity with gravity relief
            v_t = sqrt(xdot_temp(1)**2+xdot_temp(2)**2)
            a_c = v_t**2/(R_E_English - z)
            g_minus_ac = g - a_c
            A = (S_theta + (q*w-r*v)/g_minus_ac)*S_alpha
            parenthesis_term = (C_theta + cos(phi_guess)*S_theta*w/u)/&
            ((load_factor-A)/(C_alpha)+p*v/g_minus_ac)-w*S_theta/(u*C_theta)
            phi_current = acos(max(-1.0, min(1.0,  parenthesis_term)))
            current_error = abs(phi_current-phi_guess)
            i = i + 1
            phi_guess = phi_current 
        end do
        if (i >= 500) then 
            write(*,*) "Unable to converge to phi after 500 iterations given a load factor"
            write(*,*) "Current error"
            write(*,*) current_error
        end if 
        return_phi = phi_guess 
    end function calc_phi_from_load_factor

    function calc_load_factor(this, Fxb, Fzb, alpha, mass, g, ac) result(load_factor)
        implicit none 
        type(vehicle_t) :: this
        real, intent(in) :: Fxb, Fzb, alpha, mass, g, ac
        real :: load_factor
        load_factor = (Fxb*sin(alpha)-Fzb*cos(alpha))/(mass*(g-ac))
    end function calc_load_factor

    function calc_residual(this, state, p, q, r, is_trim_sideslip_angle, trim_azimuth_angle, &
        trim_bank_angle, trim_elevation_angle, trim_sideslip_angle, trim_type) result(return_state) !!! move pqr out of loop. 
        implicit none 
        type(vehicle_t) :: this
        real, intent(in) :: state(6), p, q, r, trim_azimuth_angle, trim_bank_angle, trim_elevation_angle
        real, intent(in) :: trim_sideslip_angle
        logical, intent(in) :: is_trim_sideslip_angle
        character(len=*), intent(in) :: trim_type
        real :: return_state(6)
        real :: full_state_temp(13), full_state(13)
        real :: phi, theta, psi, e0, ex, ey, ez
        real :: quaternion(4)
        real :: alpha, beta,  da, de, dr
        real :: tau, u, v, w, gravity
        real :: x,y,z,sct_pqr_coeff, temp_lhs, temp_rhs_1, temp_rhs_2
        z = this%init_state(9)
        gravity = gravity_English(-z)
        x = 0.0
        y = 0.0
        ! Diagnostic output
        if (abs(state(1)) > 1.0) then  ! alpha > ~57 degrees is suspicious
            write(*,*) "WARNING: Large alpha in calc_residual: alpha[deg]=", state(1)*180.0/PI
            write(*,*) "  Full state:", state
        end if
        alpha = state(1)
        if (trim_type == "shss" .and. is_trim_sideslip_angle) then 
            beta = trim_sideslip_angle 
            phi = state(2)
        else
            beta = state(2)
            phi = trim_bank_angle
        end if 
        ! beta = state(2)
        da = state(3)
        de = state(4)
        dr = state(5)
        tau = state(6)
        ! if (tau < 0.0) then 
            ! tau = 0.0
        ! if (tau > 1.0) then 
        !     tau = 1.0
        ! end if 
        u = this%init_V*cos(alpha)*cos(beta)
        v = this%init_V*sin(beta)
        w = this%init_V*sin(alpha)*cos(beta) 
        psi = trim_azimuth_angle
        theta = trim_elevation_angle
        this%controls(1) = da
        this%controls(2) = de 
        this%controls(3) = dr 
        this%controls(4) = tau

        quaternion = euler_to_quat([phi,theta,psi])
        e0 = quaternion(1)
        ex = quaternion(2)
        ey = quaternion(3)
        ez = quaternion(4)
        full_state_temp = [u,v,w,p,q,r,x,y,z,e0,ex,ey,ez]
        ! write(*,*) "full_state_temp", full_state_temp
        full_state = differential_equations(this, 0.0, full_state_temp)
        ! write(*,*) "full_state" 
        ! write(*,*) full_state
        return_state(1:6) = full_state(1:6)
        ! Uncomment for detailed residual debugging:
        ! write(*,'(A,6E15.6)') "  Residuals [udot,vdot,wdot,pdot,qdot,rdot]: ", return_state
    end function calc_residual

    function create_jacobian(this, states, step_size, p, q, r, is_trim_sideslip_angle, &
        trim_azimuth_angle, trim_bank_angle, trim_elevation_angle, trim_sideslip_angle, trim_type) result(jacobian) ! states should be six elements long (alpha,phi, da de dr tau) or (alpha, beta, da, de, dr, tau)
        implicit none 
        type(vehicle_t) :: this
        real, intent(in) :: step_size, p, q, r, trim_azimuth_angle, trim_bank_angle, trim_elevation_angle
        logical, intent(in) :: is_trim_sideslip_angle
        character(len=*), intent(in) :: trim_type
        real, intent(in) :: trim_sideslip_angle
        real :: states(6)
        ! real :: states_plus(6), states_minus(6)
        real :: jacobian(6,6)
        real :: R_plus(6)
        real :: R_minus(6)
        integer :: j, i
        ! write(*,'(a)') 'Building Jacobian Matrix:'
        do j = 1, 6
            states(j) = states(j) + step_size
            ! write(*,'(a,i3)') 'Computing gradient relative to G[', j-1, ']'
            ! write(*,'(a)') '   Positive Finite-Difference Step '
            ! write(*,'(a,1x,6(e25.16,","))') 'G = ', states
                ! function calc_residual(this, state, p, q, r, is_trim_sideslip_angle, trim_azimuth_angle, &
        ! trim_bank_angle, trim_elevation_angle, trim_sideslip_angle, trim_type) result(return_state) !!! move pqr out of loop. 

            R_plus = calc_residual(this, states, p, q, r, is_trim_sideslip_angle, &
                    trim_azimuth_angle, trim_bank_angle, trim_elevation_angle, &
                    trim_sideslip_angle, trim_type) ! use alpha and beta at that point to calculate u,v,w,p,q,r,phi,theta,psi which make the state vector
            states(j) = states(j) - 2*step_size 
            R_minus = calc_residual(this, states, p, q, r, is_trim_sideslip_angle, &
                    trim_azimuth_angle, trim_bank_angle, trim_elevation_angle, &
                    trim_sideslip_angle, trim_type)
            ! write(*,'(a,1x,6(e25.16,","))') 'r = ', R_plus
            ! write(*,'(a)') '   Negative Finite-Difference Step '
            ! write(*,'(a,1x,6(e25.16,","))') 'G = ', states
            ! write(*,'(a,1x,6(e25.16,","))') 'r = ', R_minus
            do i = 1, 6
                jacobian(i,j) = (R_plus(i) - R_minus(i))/(2*step_size)
            end do 

            states(j) = states(j) + step_size
        end do  
    end function create_jacobian

end module vehicle_m