module sim_m
    use adams_m
    use jsonx_m
    use linalg_mod
    use micro_time_m
    use connection_m
    implicit none
    ! Variables within sim_m
    real :: mass
    real :: I(3,3)
    real :: Ixxb, Iyyb, Izzb, Ixyb, Ixzb, Iyzb
    real :: Iinv(3,3)
    real :: h_gyro(3,3)
    real :: hdot_gyro(3)
    real :: hx, hy, hz ! non-zero terms of h_gyro
    real :: FM(6)
    real :: y_init(13)
    real :: controls(4)
    real :: controls_from_connect(4)
    real :: V_initial, alpha_initial, beta_initial
    real :: Alt_initial 
    real :: rho0,Z_temp,T_temp,P_temp,a_temp,mu_temp
    real :: trim_state(13)
    logical :: rk4_verbose, is_trim_sideslip_angle, is_use_controls, is_use_controller
    real :: trim_elevation_angle, trim_sideslip_angle, trim_bank_angle
    real :: trim_azimuth_angle, p_wind, trim_climb_angle
    real, allocatable :: rollRateControl(:), bankAngleControl(:)
    real, allocatable :: pitchRateControl(:), elevationAngleControl(:)
    real, allocatable :: yawRateControl(:), velocityControl(:)
    real :: lambda_CL, lambda_CD, lambda_Cm
    real :: alpha_0CL, alpha_sCL
    real :: alpha_0CD, alpha_sCD
    real :: alpha_0Cm, alpha_sCm, Cmmin
    logical :: include_stall

    type(connection) :: graphics, connect_controls 
    type(json_value), pointer :: j_main

    ! aero coefficients 
    real, allocatable :: aero_ref_location(:), eul0(:), angular_rates(:) 
    character(len=:),allocatable :: init_type, trim_type
    character(len=:),allocatable :: is_elevation_or_climb, is_bank_or_beta_for_shss
    real :: finite_diff_step, relax_factor, newton_tol
    real :: trim_array(9)
    real :: sref, long_ref, lat_ref
    real :: CL0, CLa, CLahat, CLqbar, CLde
    real :: CDL0, CDL1, CDL2, CDS2, CDqbar, CDaqbar, CDde, CDade, CDde2
    real :: CSb, CSpbar, CSapbar, CSrbar, CSda, CSdr
    real :: Clb, Clpbar, Clrbar, Clarbar, Clda, Cldr
    real :: Cm0, Cma, Cmqbar, Cmahat, Cmde 
    real :: Cnb, Cnpbar, Cnapbar, Cnrbar, Cnda, Cnada, Cndr 
    real :: Thrust0, Ta
    real :: weight !! might switch this out for more general stuff.
    real ::  t, u_0, v_0, w_0, p_0, q_0, r_0, x_0, y_0, z_0
    real :: dt, tf, delta_t_over_2, delta_t_over_6
    real :: test_quat(4)
    real :: test_quat_2(4)
    real :: quat_mult_result(4)
    integer :: n, newton_max_iter
    contains
    ! end subroutine simulation_main
    subroutine run()
        implicit none 
        real :: y(13), y1(13), s(14)
        ! real :: Z_temp, P_temp, T_temp, a_temp, mu_temp
        real :: cpu_start_time, cpu_end_time, time1, time2, actual_time, integrated_time
        integer :: io_unit
        logical :: real_time
        delta_t_over_2 = dt/2.0
        delta_t_over_6 = dt/6.0
        n = 13
        ! call std_atm_English(0.0,Z_temp,T_temp,P_temp,rho0,a_temp,mu_temp)
        ! initial conditions 
        ! test_euler = [50.0*PI/180.0, 0.0, 0.0]
        test_quat = [1.0, 2.0, 3.0, 4.0]
        test_quat_2 = [1.0,-2.0,-3.0,-4.0]
        quat_mult_result = quat_mult(test_quat, test_quat_2)
        write(*,*) "QUAT MULT RESULT", quat_mult_result
        write(*,*) "QUAT MULT RESULT", quat_mult_result
        write(*,*) "QUAT MULT RESULT", quat_mult_result

        t = 0.0
        y = y_init
        real_time = .false.
        if (abs(dt)<TOLERANCE) then 
            real_time = .true.
            time1 = get_time()
            y1 = runge_kutta(t, y,dt)
            call quat_norm(y1(10:13))
            time2 = get_time()
            dt = time2-time1 
            write(*,*) "Dt = ", dt
            if (dt == 0.0) then 
                write(*,*) "THE DT IS EXACTLY ZERO!!!"
            end if 
            y = y_init 
            t = 0.0
        end if
        delta_t_over_2 = dt/2.0
        delta_t_over_6 = dt/6.0
        open(newunit=io_unit, file='output.txt', status='replace', action='write')
        call quat_norm(y(10:13))
        write(io_unit,*) "        t[s]                   u[ft/s]               v[ft/s]"// &
                   "                 w[ft/s]                p[rad/s]             q[rad/s]"  // &
                   "            r[rad/s]                x[ft]                  y[ft]" // &                   
                   "                  z[ft]                    e0                   ex" // &
                   "                      ey                   ez"
        write(io_unit,'(14E22.13)') t,y(:) ! 14E22.13 is for 14 numbers in scientific notation, each 22 characters wide with 13 after the decimal
        ! write(*,'(14E22.13)') t,y(:) ! 14E22.13 is for 14 numbers in scientific notation, each 22 characters wide with 13 after the decimal
        ! do while(y(9)<0.0 .and. t<tf) ! while altitude is greater than 0 ft (altitude is positive going down in our coordinate systems) or time is less than final time
        cpu_start_time = get_time()
        time1 = cpu_start_time
        integrated_time = 0.0

        do while(t<tf) ! while altitude is greater than 0 ft (altitude is positive going down in our coordinate systems) or time is less than final time
            ! write(*,*) "dt", dt
            if (is_use_controls) then 
                controls = connect_controls%recv()
            end if
            y1 = runge_kutta(t,y,dt)
            call quat_norm(y1(10:13))
            y = y1
            t = t + dt 
            integrated_time = integrated_time + dt
            ! write(io_unit,'(14E22.13)') t,y(:)
            s(1) = t
            s(2:14) = y(1:13)
            ! write(*,'(14E22.13)') s
            call graphics%send(s)
            ! write(*,'(14E22.13)') t,y(:)
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

    function cross_product_3D(a, b) result(result)
        implicit none
        real, intent(in) :: a(3), b(3) 
        real :: result(3)
        result(1) = a(2)*b(3) - a(3)*b(2)
        result(2) = a(3)*b(1) - a(1)*b(3)
        result(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product_3D

    function runge_kutta(t_0, state, delta_t) result(state_out)
        implicit none 
        real, intent(in) :: t_0
        real, intent(in), dimension(:) :: state
        real, intent(in) :: delta_t
        real :: state_out(size(state))
        real :: k1(n), k2(n), k3(n), k4(n)
        real :: state_temp(n)
        real :: t_0_plus_delta_t_over_2
        t_0_plus_delta_t_over_2 = t_0 + delta_t_over_2
        k1 = differential_equations(t_0, state)
            state_temp = state
            state_temp = state_temp + delta_t_over_2 * k1
        k2 = differential_equations(t_0_plus_delta_t_over_2, state_temp)
            state_temp = state
            state_temp = state_temp + delta_t_over_2 * k2
        k3 = differential_equations(t_0_plus_delta_t_over_2, state_temp)
            state_temp = state 
            state_temp = state_temp + delta_t * k3 
        k4 = differential_equations(t_0 + delta_t, state_temp)
        state_out = state + delta_t_over_6 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    end function runge_kutta

    function differential_equations(time, state) result(res)
        implicit none 
        real, intent(in) :: time
        real, intent(in), dimension(:) :: state
        real :: res(size(state))
        real :: u,v,w,p,q,r,x,y,z
        real :: e0,ex,ey,ez
        real, dimension(3) :: pqr_temp, rot_and_inertia_temp
        real :: gravity
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
        call pseudo_aero(state)
        rot_and_inertia_temp = & 
        [FM(4) + dot_product(h_gyro(1,:),pqr_temp) + ((Iyyb-Izzb)*q*r-Iyzb*(q**2-r**2)-Ixzb*p*q+Ixyb*p*r)-hdot_gyro(1),&
        FM(5) + dot_product(h_gyro(2,:),pqr_temp) + ((Izzb-Ixxb)*p*r-Ixzb*(r**2-p**2)-Ixyb*q*r+Iyzb*p*q)-hdot_gyro(2),& 
        FM(6) + dot_product(h_gyro(3,:),pqr_temp) + ((Ixxb-Iyyb)*p*q-Ixyb*(p**2-q**2)-Iyzb*p*r+Ixzb*q*r)-hdot_gyro(3)]
        res(1) = 1/mass * FM(1) + gravity * 2 *(ex*ez-ey*e0) + r*v - q*w ! udot body-fixed
        res(2) = 1/mass * FM(2) + gravity * 2 *(ey*ez+ex*e0) + p*w - r*u ! vdot body-fixed
        res(3) = 1/mass * FM(3) + gravity * (ez**2+e0**2-ex**2-ey**2) + q*u - p*v ! wdot body-fixed
        res(4:6) = matmul(Iinv, rot_and_inertia_temp) ! pdot, qdot, rdot body-fixed
        res(7:9) = quat_dependent_to_base((/u,v,w/), (/e0, ex, ey, ez/)) !!! + wind ! xdot ydot zdot earth-fixed
        res(10) = 0.5 * dot_product((/-ex, -ey, -ez/),pqr_temp) !e0
        res(11) = 0.5 * dot_product((/e0, -ez, ey/),pqr_temp) !ex
        res(12) = 0.5 * dot_product((/ez, e0, -ex/),pqr_temp) !ey
        res(13) = 0.5 * dot_product((/-ey, ex, e0/),pqr_temp) !ez
        if (rk4_verbose) then
            write(*,*) "state"
            write(*,*) state
            write(*,*) ""
            write(*,*) "res"
            write(*,'(14E22.13)') res
            write(*,*) ""
        end if
    end function differential_equations

    subroutine init(filename)
        implicit none 
        character(100), intent(in) :: filename
        type(json_value), pointer :: j_connections, j_graphics, j_controls
        ! type2, intent(out) ::  arg2
        ! call get_command_argument(1,filename)
        call std_atm_English(0.0,Z_temp,T_temp,P_temp,rho0,a_temp,mu_temp)

        call jsonx_load(filename,j_main)
        ! simulation
        call jsonx_get(j_main, "simulation.time_step[sec]", dt, 0.0)
        call jsonx_get(j_main, "simulation.end_time[sec]", tf, 20.0)
        call jsonx_get(j_main, "simulation.rk4_verbose", rk4_verbose, .false.)
        ! vehicle
        ! thrust
        call jsonx_get(j_main, "vehicle.thrust.Thrust0[lbf]", Thrust0)
        call jsonx_get(j_main, "vehicle.thrust.Ta", Ta)
        ! aerodynamics
        ! reference
        call jsonx_get(j_main, "vehicle.aerodynamics.reference.area[ft^2]", sref)
        call jsonx_get(j_main, "vehicle.aerodynamics.reference.longitudinal_length[ft]", long_ref)
        call jsonx_get(j_main, "vehicle.aerodynamics.reference.lateral_length[ft]", lat_ref)
        call jsonx_get(j_main, "vehicle.aerodynamics.reference.relative_location[ft]", aero_ref_location,0.0,3)
        ! coefficients 
        !CL 
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CL.0", CL0)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CL.alpha", CLa)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CL.alphahat", CLahat)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CL.qbar", CLqbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CL.elevator", CLde)
        !CS 
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.beta", CSb)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.pbar", CSpbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.alpha_pbar", CSapbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.rbar", CSrbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.aileron", CSda)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CS.rudder", CSdr)
        ! CD
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.L0", CDL0)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.CL1", CDL1)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.CL1_CL1", CDL2)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.CS_CS", CDS2)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.qbar", CDqbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.alpha_qbar", CDaqbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.elevator", CDde)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.alpha_elevator", CDade)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.CD.elevator_elevator", CDde2)
        ! Cl
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.beta", Clb)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.pbar", Clpbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.rbar", Clrbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.alpha_rbar", Clarbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.aileron", Clda)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cl.rudder", Cldr)
        ! Cm 
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cm.0", Cm0)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cm.alpha", Cma)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cm.qbar", Cmqbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cm.alphahat", Cmahat)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cm.elevator", Cmde)
        ! Cn
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.beta", Cnb)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.pbar", Cnpbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.alpha_pbar", Cnapbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.rbar", Cnrbar)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.aileron", Cnda)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.alpha_aileron", Cnada)
        call jsonx_get(j_main, "vehicle.aerodynamics.coefficients.Cn.rudder", Cndr)
        ! initial conditions
        call jsonx_get(j_main, "connections.controls.is_use", is_use_controls)
        call mass_inertia() ! computes mass, inertia, and gyroscopic properties of a thing which is weighed at a geopotential altitude of exactly zero. 
        y_init = 0.0
        call jsonx_get(j_main, "initial.time[sec]", t, 0.0)
        call jsonx_get(j_main, "initial.airspeed[ft/sec]", V_initial)
        call jsonx_get(j_main, "initial.altitude[ft]", y_init(9))
        y_init(9) = - y_init(9)
        Alt_initial = y_init(9)
        call jsonx_get(j_main, "initial.Euler_angles[deg]", eul0,0.0,3)
        eul0 = eul0*PI/180.0
        y_init(10:13) = euler_to_quat(eul0)
        call jsonx_get(j_main, "initial.type", init_type)
        is_bank_or_beta_for_shss = "bank"
        if (init_type=="state") then
            call jsonx_get(j_main, "initial.state.angle_of_attack[deg]", alpha_initial)
            alpha_initial = alpha_initial * PI/180.0
            call jsonx_get(j_main, "initial.state.sideslip_angle[deg]", beta_initial)
            beta_initial = beta_initial * PI/180.0
            y_init(1) = V_initial*cos(alpha_initial)*cos(beta_initial)
            y_init(2) = V_initial*sin(beta_initial)
            y_init(3) = V_initial*sin(alpha_initial)*cos(beta_initial)
            call jsonx_get(j_main, "initial.state.angular_rates[deg/s]", angular_rates,0.0,3)
            y_init(4:6) = angular_rates * PI/180.0
            call jsonx_get(j_main, "initial.state.aileron[deg]", controls(1))
            call jsonx_get(j_main, "initial.state.elevator[deg]", controls(2))
            call jsonx_get(j_main, "initial.state.rudder[deg]", controls(3))
            call jsonx_get(j_main, "initial.state.throttle", controls(4))
            controls(1:3) = controls(1:3)*PI/180.0
        else 
            call jsonx_get(j_main, "initial.trim.type", trim_type)
            call jsonx_get(j_main, "initial.trim.elevation_or_climb", is_elevation_or_climb)
            if (is_elevation_or_climb == "climb") then
                call jsonx_get(j_main, "initial.trim.climb_angle[deg]", trim_climb_angle)
                trim_climb_angle = trim_climb_angle*PI/180.0
                trim_elevation_angle = 0.0
            else 
                call jsonx_get(j_main, "initial.trim.elevation_angle[deg]", trim_elevation_angle)
                trim_elevation_angle = trim_elevation_angle*PI/180.0
                trim_climb_angle = 0.0
            end if 
            call jsonx_get(j_main, "initial.trim.solver.finite_difference_step_size", finite_diff_step)
            call jsonx_get(j_main, "initial.trim.solver.relaxation_factor", relax_factor)
            call jsonx_get(j_main, "initial.trim.solver.tolerance", newton_tol)
            call jsonx_get(j_main, "initial.trim.solver.max_iterations", newton_max_iter)
            trim_bank_angle = 0.0
            trim_sideslip_angle = 0.0
            p_wind = 0.0
            trim_azimuth_angle = eul0(3)
            is_trim_sideslip_angle = .false.
            if (trim_type == "sct") then
                call jsonx_get(j_main, "initial.trim.type_sct.bank_angle[deg]", trim_bank_angle)
            else if (trim_type == "shss") then
                call jsonx_get(j_main, "initial.trim.type_shss.bank_or_beta", is_bank_or_beta_for_shss)
                if (is_bank_or_beta_for_shss == "bank") then 
                    call jsonx_get(j_main, "initial.trim.type_shss.bank.bank_angle[deg]", trim_bank_angle)
                else if (is_bank_or_beta_for_shss == "beta") then 
                    is_trim_sideslip_angle = .true.
                    call jsonx_get(j_main, "initial.trim.type_shss.beta.sideslip_angle[deg]", trim_sideslip_angle)
                end if 
            else if (trim_type == "vbr") then 
                call jsonx_get(j_main, "initial.trim.type_vbr.p_wind[deg/s]", p_wind)
                call jsonx_get(j_main, "initial.trim.type_vbr.bank_angle[deg]", trim_bank_angle)
            end if 
            trim_bank_angle = trim_bank_angle*PI/180.0
            trim_sideslip_angle = trim_sideslip_angle*PI/180.0
            p_wind = p_wind*PI/180.0
            ! call jsonx_get(j_main, "initial.trim.elevation_angle[deg]", trim_elevation_angle)
            ! call jsonx_get(j_main, "initial.trim.bank_angle[deg]", trim_bank_angle)
            ! call jsonx_get(j_main, "initial.trim.sideslip_angle[deg]", trim_sideslip_angle)
            alpha_initial = 0.0
            beta_initial = 0.0
            if (trim_type == "shss" .and. is_trim_sideslip_angle) then
                beta_initial = trim_sideslip_angle
            end if 
            write(*,*) 
            write(*,*) 
            write(*,*) "alpha", alpha_initial
            write(*,*) "beta", beta_initial
            write(*,*) "x[ft]", y_init(7)
            write(*,*) "y[ft]", y_init(8)
            write(*,*) "altitude[ft]", y_init(9)
            write(*,*) "trim_type ", trim_type
            write(*,*) "is_bank_or_beta_for_shss ", is_bank_or_beta_for_shss
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
            write(*,'(a,1x,e25.16)') 'Initial beta [deg]  = ', beta_initial*180.0/PI

            write(*,'(a)') 'Newton Solver Settings:'
            write(*,'(a,1x,e25.16)') 'Finite Difference Step Size = ', finite_diff_step
            write(*,'(a,1x,e25.16)') '          Relaxation Factor = ', relax_factor
            write(*,'(a,1x,e25.16)') '                  Tolerance = ', newton_tol

            trim_array = trim_algorithm(y_init(9), newton_tol)
            if (is_bank_or_beta_for_shss == "beta" .and. trim_type == "shss") then 
                y_init(10:13) = euler_to_quat([trim_array(2), trim_elevation_angle, trim_azimuth_angle])
                write(*,'(A12,1X,E22.13)') "theta[deg]", trim_elevation_angle*180.0/PI
                write(*,'(A12,1X,E22.13)') "phi[deg]", trim_array(2)*180.0/PI
                write(*,'(A12,1X,E22.13)') "alpha[deg]", trim_array(1)*180.0/PI
                write(*,'(A12,1X,E22.13)') "beta[deg]", beta_initial*180.0/PI
                write(*,'(A12,1X,E22.13)') "p[deg/s]", trim_array(3)*180.0/PI
                write(*,'(A12,1X,E22.13)') "q[deg/s]", trim_array(4)*180.0/PI
                write(*,'(A12,1X,E22.13)') "r[deg/s]", trim_array(5)*180.0/PI
                write(*,'(A12,1X,E22.13)') "da[deg]", trim_array(6)*180.0/PI
                write(*,'(A12,1X,E22.13)') "de[deg]", trim_array(7)*180.0/PI
                write(*,'(A12,1X,E22.13)') "dr[deg]", trim_array(8)*180.0/PI
                write(*,'(A12,1X,E22.13)') "tau", trim_array(9)
            else
                beta_initial = trim_array(2)
                y_init(10:13) = euler_to_quat([trim_bank_angle, trim_elevation_angle, trim_azimuth_angle])
                write(*,'(A12,1X,E22.13)') "theta[deg]", trim_elevation_angle*180.0/PI
                write(*,'(A12,1X,E22.13)') "phi[deg]", trim_bank_angle*180.0/PI
                write(*,'(A12,1X,E22.13)') "alpha[deg]", trim_array(1)*180.0/PI
                write(*,'(A12,1X,E22.13)') "beta[deg]", trim_array(2)*180.0/PI
                write(*,'(A12,1X,E22.13)') "p[deg/s]", trim_array(3)*180.0/PI
                write(*,'(A12,1X,E22.13)') "q[deg/s]", trim_array(4)*180.0/PI
                write(*,'(A12,1X,E22.13)') "r[deg/s]", trim_array(5)*180.0/PI
                write(*,'(A12,1X,E22.13)') "da[deg]", trim_array(6)*180.0/PI
                write(*,'(A12,1X,E22.13)') "de[deg]", trim_array(7)*180.0/PI
                write(*,'(A12,1X,E22.13)') "dr[deg]", trim_array(8)*180.0/PI
                write(*,'(A12,1X,E22.13)') "tau", trim_array(9)
                ! write(*,'(A12,1X,E22.13)') "psi[deg]", trim_azimuth_angle*180.0/PI
            end if  
            alpha_initial = trim_array(1)
            y_init(1) = V_initial*cos(alpha_initial)*cos(beta_initial)
            y_init(2) = V_initial*sin(beta_initial)
            y_init(3) = V_initial*sin(alpha_initial)*cos(beta_initial)
            y_init(4:6) = trim_array(3:5)
            controls(1:4) = trim_array(6:9) 
        end if
        ! Stall Stuff 
        call jsonx_get(j_main, "vehicle.aerodynamics.stall.include_stall", include_stall)
        call jsonx_get(j_main, "vehicle.aerodynamics.stall.CL.alpha_0[deg]", alpha_0CL)
        call jsonx_get(j_main, "vehicle.aerodynamics.stall.CL.alpha_s[deg]", alpha_sCL)
        call jsonx_get(j_main, "vehicle.aerodynamics.stall.CL.lambda_b", lambda_CL)
        write(*,*) "Is Include Stall? ", include_stall
        write(*,*) "alpha_0CL", alpha_0CL
        write(*,*) "alpha_sCL", alpha_sCL
        write(*,*) "lambda_CL", lambda_CL
        alpha_0CL = alpha_0CL*PI/180.0
        alpha_sCL = alpha_sCL*PI/180.0
        call jsonx_get(j_main, "vehicle.aerodynamics.stall.CD.alpha_0[deg]", alpha_0CD)
        call jsonx_get(j_main, "vehicle.aerodynamics.stall.CD.alpha_s[deg]", alpha_sCD)
        call jsonx_get(j_main, "vehicle.aerodynamics.stall.CD.lambda_b", lambda_CD)
        write(*,*) "alpha_0CD", alpha_0CD
        write(*,*) "alpha_sCD", alpha_sCD
        write(*,*) "lambda_CD", lambda_CD
        alpha_0CD = alpha_0CD*PI/180.0
        alpha_sCD = alpha_sCD*PI/180.0

        call jsonx_get(j_main, "vehicle.aerodynamics.stall.Cm.min", Cmmin)
        call jsonx_get(j_main, "vehicle.aerodynamics.stall.Cm.alpha_0[deg]", alpha_0Cm)
        call jsonx_get(j_main, "vehicle.aerodynamics.stall.Cm.alpha_s[deg]", alpha_sCm)
        call jsonx_get(j_main, "vehicle.aerodynamics.stall.Cm.lambda_b", lambda_Cm)
        write(*,*) "Cm min ", Cmmin
        write(*,*) "alpha_0Cm", alpha_0Cm
        write(*,*) "alpha_sCm", alpha_sCm
        write(*,*) "lambda_Cm", lambda_Cm
        alpha_0Cm = alpha_0Cm*PI/180.0
        alpha_sCm = alpha_sCm*PI/180.0
        ! controller stuff 
        call jsonx_get(j_main, "controller.is_use_controller", is_use_controller)
        call jsonx_get(j_main, "controller.rollRateControl", rollRateControl,0.0,4)
        call jsonx_get(j_main, "controller.bankAngleControl", bankAngleControl,0.0,4) 
        call jsonx_get(j_main, "controller.pitchRateControl", pitchRateControl,0.0,4) 
        call jsonx_get(j_main, "controller.elevationAngleControl", elevationAngleControl,0.0,4)
        call jsonx_get(j_main, "controller.yawRateControl", yawRateControl,0.0,4) 
        call jsonx_get(j_main, "controller.velocityControl", velocityControl,0.0,4)
        ! connections 
        call jsonx_get(j_main, 'connections', j_connections)
        call jsonx_get(j_connections, 'graphics', j_graphics)
        call graphics%init(j_graphics)
        call jsonx_get(j_connections, 'controls', j_controls)
        call connect_controls%init(j_controls)
        call print_aero_table()
    end subroutine init

    function trim_algorithm(H_altitude, newton_tol) result(trim_result)
        implicit none
        real, intent(in) :: H_altitude, newton_tol
        real :: trim_result(9)
        real :: alpha, beta, p, q, r, da, de, dr, tau,x,y,z
        real :: pos(3), quat_orientation(4)
        real :: phi, theta, psi
        real :: current_error
        real :: u,v,w, sct_pqr_coeff
        real, allocatable :: DeltaG(:)
        real :: residual(6)
        real :: jacobian(6,6)
        real :: gravity 
        real :: newton_input(6)
        integer :: i, j
        gravity = gravity_English(-H_altitude)
        alpha = 0.0
        allocate(DeltaG(6))
        if (trim_type == "shss" .and. is_trim_sideslip_angle) then 
            beta = trim_sideslip_angle
        else
            beta = 0.0
        end if 
        u = V_initial*cos(alpha)*cos(beta)
        v = V_initial*sin(beta)
        w = V_initial*sin(alpha)*cos(beta) 
        x = 0.0
        y = 0.0
        z = H_altitude
        da = 0.0
        de = 0.0 
        dr = 0.0
        tau = 0.0
        if (trim_type == "shss") then
            if (is_trim_sideslip_angle) then 
                beta = trim_sideslip_angle
            else 
                phi = trim_bank_angle
            end if 
        else
            phi = trim_bank_angle
        end if 
        ! trim_elevation_angle = calc_theta_from_climb_angle(trim_climb_angle, u, v, w, phi)
        if (trim_type == "sct") then 
            sct_pqr_coeff = gravity*sin(trim_bank_angle)*cos(trim_elevation_angle)/&
            (u*cos(trim_elevation_angle)*cos(trim_bank_angle)+w*sin(trim_elevation_angle))
            p = -sct_pqr_coeff*(sin(trim_elevation_angle))
            q = sct_pqr_coeff*(sin(trim_bank_angle)*cos(trim_elevation_angle))
            r = sct_pqr_coeff*(cos(trim_bank_angle)*cos(trim_elevation_angle))
        else if (trim_type == "vbr") then 
            p = (p_wind/V_initial)*u
            q = (p_wind/V_initial)*v
            r = (p_wind/V_initial)*w
        else
            p = 0.0
            q = 0.0
            r = 0.0
        end if 
        current_error = 100.0
        if (trim_type == "shss" .and. is_trim_sideslip_angle) then
                newton_input = [alpha, phi, da, de, dr, tau]
            else
                newton_input = [alpha, beta, da, de, dr, tau]
        end if
        j = 1
        do while(current_error > newton_tol)
            if (newton_input(6) < 0.0) then 
                newton_input(6) = 0.0
            else if (newton_input(6) > 1.0) then 
                newton_input(6) = 1.0
            end if 
            alpha = newton_input(1)
            if (trim_type == "shss" .and. is_trim_sideslip_angle) then 
                beta = trim_sideslip_angle
                phi = newton_input(2)
            else
                beta = newton_input(2)
            end if 
            ! beta = newton_input(2)
            u = V_initial*cos(alpha)*cos(beta)
            v = V_initial*sin(beta)
            w = V_initial*sin(alpha)*cos(beta) 
            if (is_elevation_or_climb == "climb") then 
                trim_elevation_angle = calc_theta_from_climb_angle(trim_climb_angle, u, v, w, phi)
            end if 
            ! write(*,*) "Trim theta", trim_elevation_angle * 180.0/PI
            if (trim_type == "sct") then  
                sct_pqr_coeff = gravity*sin(trim_bank_angle)*cos(trim_elevation_angle)/&
                (u*cos(trim_elevation_angle)*cos(trim_bank_angle)+w*sin(trim_elevation_angle))
                p = -sct_pqr_coeff*(sin(trim_elevation_angle))
                q = sct_pqr_coeff*(sin(trim_bank_angle)*cos(trim_elevation_angle))
                r = sct_pqr_coeff*(cos(trim_bank_angle)*cos(trim_elevation_angle))
            else if (trim_type == "vbr") then 
                p = (p_wind/V_initial)*u
                q = (p_wind/V_initial)*v
                r = (p_wind/V_initial)*w
            else 
                p = 0.0
                q = 0.0
                r = 0.0
            end if 
            residual = calc_residual(newton_input, p, q, r)
            ! if (rk4_verbose) then
            ! write(*,'(a)') 'Updating rotation rates for ', trim_type
            ! write(*,'(a,1x,e25.16)') 'p [deg/s] = ', (p * 180.0/PI)
            ! write(*,'(a,1x,e25.16)') 'q [deg/s] = ', (q * 180.0/PI)
            ! write(*,'(a,1x,e25.16)') 'r [deg/s] = ', (r * 180.0/PI)
            ! write(*,'(a)') 'G defined as G = [alpha, beta, aileron, elevator, rudder, throttle]'
            ! write(*,'(a,1x,6(e25.16,","))') 'G = ', newton_input
            ! write(*,'(a,1x,6(e25.16,","))') 'r = ', residual
            ! write(*,'(a,1x,e25.16)') 'current_error = ', current_error
            ! end if
            ! write(*,*) "residual"
            ! write(*,*) residual
            jacobian = create_jacobian(newton_input, finite_diff_step, p, q, r)
            ! if (rk4_verbose) then
            ! write(*,'(a)') 'Jacobian Matrix ='
            ! do i = 1,6
            !     write(*,'(6(1x,e25.16))') jacobian(i,:)
            ! end do
            ! end if
            ! write(*,*) "jacobian"
            ! write(*,*) jacobian
            call lu_solve(6,jacobian,residual,DeltaG)
            ! write(*,*) "DeltaG"
            ! write(*,*) -DeltaG
            newton_input = newton_input - relax_factor*DeltaG
            ! write(*,*) "New G"
            ! write(*,*) newton_input
            residual = calc_residual(newton_input, p, q, r)
            ! write(*,'(a,1x,6(e25.16,","))') 'r = ', residual
            current_error = maxval(abs(residual))
            ! write(*,*) "Iteration, Residual, alpha, beta, p, q, &
            ! r, phi, theta, ail, el, rud, throttle"
            ! write(*,*) j, current_error, newton_input(1)*180.0/PI, &
            ! newton_input(2)*180.0/PI, p*180.0/PI, q*180.0/PI, r*180.0/PI, phi*180.0/PI,trim_elevation_angle*180.0/PI, &
            ! newton_input(3)*180.0/PI, newton_input(4)*180.0/PI, newton_input(5)*180.0/PI, newton_input(6)
            j = j + 1
        end do 
        alpha = newton_input(1)
        beta = newton_input(2)
        da = newton_input(3)
        de = newton_input(4)
        dr = newton_input(5)
        tau = newton_input(6)
        u = V_initial*cos(alpha)*cos(beta)
        v = V_initial*sin(beta)
        w = V_initial*sin(alpha)*cos(beta) 
        if (trim_type == "shss" .and. is_trim_sideslip_angle) then
            phi = newton_input(2)
        else 
            phi = trim_bank_angle
        end if
        theta = trim_elevation_angle
        psi = trim_azimuth_angle
        if (trim_type == "sct") then 
            sct_pqr_coeff = gravity*sin(phi)*cos(theta)/(u*cos(theta)*cos(phi)+w*sin(theta))
            p = -sct_pqr_coeff*(sin(theta))
            q = sct_pqr_coeff*(sin(phi)*cos(theta))
            r = sct_pqr_coeff*(cos(phi)*cos(theta))
        else if (trim_type == "vbr") then 
            p = (p_wind/V_initial)*u
            q = (p_wind/V_initial)*v
            r = (p_wind/V_initial)*w
        else 
            p = 0.0
            q = 0.0
            r = 0.0
        end if 
        trim_result = [alpha, beta, p, q, r, da, de, dr, tau]
    end function trim_algorithm

    function create_jacobian(states, step_size, p, q, r) result(jacobian) ! states should be six elements long (alpha,phi, da de dr tau) or (alpha, beta, da, de, dr, tau)
        implicit none 
        real, intent(in) :: step_size, p, q, r
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
            R_plus = calc_residual(states, p, q, r) ! use alpha and beta at that point to calculate u,v,w,p,q,r,phi,theta,psi which make the state vector
            states(j) = states(j) - 2*step_size 
            R_minus = calc_residual(states, p, q, r)
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

    function calc_theta_from_climb_angle(climb_angle, u, v, w, phi) result(return_theta)
        implicit none 
        real, intent(in) :: climb_angle, u, v, w, phi
        real :: return_theta, error_for_elevation, temp_lhs, temp_rhs_1, temp_rhs_2
        real :: parenth_temp, S_theta_1, S_theta_2, theta_1, theta_2
        ! write(*,*) "Solving for elevation angle given a climb angle"
        error_for_elevation = 1e-12
        temp_lhs = V_initial*sin(trim_climb_angle)
        parenth_temp = v*sin(phi)+w*cos(phi)
        S_theta_1 = (u*V_initial*sin(trim_climb_angle)+parenth_temp*sqrt(u**2 + parenth_temp**2 - &
        V_initial**2*sin(trim_climb_angle)**2))/(u**2 + parenth_temp**2)
        theta_1 = asin(S_theta_1)
        ! write(*,*) "theta 1 [deg] = ", theta_1 * 180.0/PI
        S_theta_2 = (u*V_initial*sin(trim_climb_angle)-parenth_temp*sqrt(u**2 + parenth_temp**2 - &
        V_initial**2*sin(trim_climb_angle)**2))/(u**2 + parenth_temp**2)
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
            write(*,*) "WARNING, BOTH THETA VALUES DO NOT SATISFY THE LHS OF EQ. 7.2.3" 
        end if 
    end function calc_theta_from_climb_angle

    function calc_residual(state, p, q, r) result(return_state) !!! move pqr out of loop. 
        implicit none 
        real, intent(in) :: state(6), p, q, r
        real :: return_state(6)
        real :: full_state_temp(13), full_state(13)
        real :: phi, theta, psi, e0, ex, ey, ez
        real :: quaternion(4)
        real :: alpha, beta,  da, de, dr
        real :: tau, u, v, w, gravity
        real :: x,y,z,sct_pqr_coeff, temp_lhs, temp_rhs_1, temp_rhs_2
        z = y_init(9)
        gravity = gravity_English(-z)
        x = 0.0
        y = 0.0
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
        if (tau < 0.0) then 
            tau = 0.0
        else if (tau > 1.0) then 
            tau = 1.0
        end if 
        u = V_initial*cos(alpha)*cos(beta)
        v = V_initial*sin(beta)
        w = V_initial*sin(alpha)*cos(beta) 
        psi = trim_azimuth_angle
        theta = trim_elevation_angle
        controls(1) = da
        controls(2) = de 
        controls(3) = dr 
        controls(4) = tau

        ! if (trim_type == "shss") then
        !     if (is_trim_sideslip_angle) then 
        !         beta = trim_sideslip_angle
        !         phi = state(2)
        !         u = V_initial*cos(alpha)*cos(beta)
        !         v = V_initial*sin(beta)
        !         w = V_initial*sin(alpha)*cos(beta) 
        !     else 
        !         phi = trim_bank_angle
        !     end if 
        ! end if 
        quaternion = euler_to_quat([phi,theta,psi])
        e0 = quaternion(1)
        ex = quaternion(2)
        ey = quaternion(3)
        ez = quaternion(4)
        full_state_temp = [u,v,w,p,q,r,x,y,z,e0,ex,ey,ez]
        ! write(*,*) "full_state_temp", full_state_temp
        full_state = differential_equations(0.0, full_state_temp)
        ! write(*,*) "full_state" 
        ! write(*,*) full_state
        return_state(1:6) = full_state(1:6)
    end function calc_residual

    subroutine mass_inertia() 
        implicit none 
        real :: gravity
        real :: I(3,3)
        real :: det_I
        real :: I_tilde(3,3) 
        real :: a11, a22, a33, a12, a13, a21, a23, a31, a32
        I = 0.0
        call jsonx_get(j_main, "vehicle.mass.weight[lbf]", weight)
        call jsonx_get(j_main, "vehicle.mass.Ixx[slug-ft^2]", I(1,1))
        call jsonx_get(j_main, "vehicle.mass.Iyy[slug-ft^2]", I(2,2))
        call jsonx_get(j_main, "vehicle.mass.Izz[slug-ft^2]", I(3,3))
        call jsonx_get(j_main, "vehicle.mass.Ixy[slug-ft^2]", I(1,2))
        call jsonx_get(j_main, "vehicle.mass.Ixz[slug-ft^2]", I(3,1))
        call jsonx_get(j_main, "vehicle.mass.Iyz[slug-ft^2]", I(3,2))
        call jsonx_get(j_main, "vehicle.mass.hx[slug-ft^2/s]", hx) ! in the 3 by 3 h matrix in Eq. 5.4.6
        call jsonx_get(j_main, "vehicle.mass.hy[slug-ft^2/s]", hy)
        call jsonx_get(j_main, "vehicle.mass.hz[slug-ft^2/s]", hz)
        I(1,2) = -I(1,2)
        I(3,1) = -I(3,1)
        I(3,2) = -I(3,2)
        I(2,1) = I(1,2)
        I(1,3) = I(3,1)
        I(2,3) = I(3,2)
        Ixxb = I(1,1)
        Iyyb = I(2,2)
        Izzb = I(3,3)
        Ixyb = I(1,2)
        Ixzb = I(1,3)
        Iyzb = I(2,3)
        a11 = I(1,1)
        a22 = I(2,2)
        a33 = I(3,3)
        a12 = I(1,2)
        a13 = I(1,3)
        a21 = I(2,1)
        a23 = I(2,3)
        a31 = I(3,1)
        a32 = I(3,2)
        gravity = gravity_English(0.0)
        mass = weight/gravity
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
        Iinv = 1/(det_I)*I_tilde
        h_gyro = 0.0
        h_gyro(1,2) = -hz
        h_gyro(1,3) = hy
        h_gyro(2,1) = hz
        h_gyro(2,3) = -hx
        h_gyro(3,1) = -hy
        h_gyro(3,2) = hx
        hdot_gyro = 0.0
    end subroutine mass_inertia

    subroutine pseudo_aero(y) 
        implicit none
        real, intent(in) :: y(13)
        real :: da, de, dr, tau
        real :: V, alpha, beta, pbar, qbar, rbar, ahat
        real :: CL1, CL, CS, CD, Cll, Cm, Cn 
        real :: sa, ca, sb, cb, sign_a
        real :: Z, T, P, rho, a, mu
        real :: exp_pos, exp_neg, sigma
        real :: CL_newt, CD_newt, Cm_newt
        real :: CL1_blended, CD1_blended, Cm1_blended
        ahat = 0.0
        !!!! receive controls from python script here !!!!
        da = controls(1)
        de = controls(2)
        dr = controls(3)
        tau = controls(4)
        if (tau < 0.0) then 
            tau = 0.0
        else if (tau > 1.0) then 
            tau = 1.0
        end if 

        call std_atm_English(-y(9), Z, T, P, rho, a, mu)

        V = sqrt(y(1)**2 + y(2)**2 + y(3)**2)
        alpha = atan2(y(3), y(1)) ! Eq. 3.4.4
        beta = asin(y(2)/V) ! Eq. 3.4.5
        pbar = 0.5*y(4)*lat_ref/(V)
        qbar = 0.5*y(5)*long_ref/(V)
        rbar = 0.5*y(6)*lat_ref/(V)

        sa = sin(alpha)
        ca = cos(alpha)
        sb = sin(beta)
        cb = cos(beta)
        sign_a = sign(1.0, alpha)

        CL1 = CL0 +CLa*alpha
        CL = CL1 + CLqbar*qbar+CLahat*ahat + CLde*de
        CS = CSb*beta + (CSpbar+CSapbar*alpha)*pbar + CSrbar*rbar + CSda*da + CSdr*dr
        CD = CDL0 + CDL1*CL1 + CDL2*CL1**2 + CDS2*CS**2 + (CDqbar + CDaqbar*alpha)*qbar + (CDde + CDade*alpha)*de + CDde2*de**2
        Cll = Clb*beta + Clpbar*pbar + (Clrbar + Clarbar*alpha)*rbar + Clda*da + Cldr*dr
        Cm = Cm0 + Cma*alpha + Cmqbar*qbar + Cmahat*ahat + Cmde*de 
        Cn = Cnb*beta + (Cnpbar + Cnapbar*alpha)*pbar + Cnrbar*rbar + (Cnda + Cnada*alpha)*da + Cndr*dr
        
        if (include_stall) then
            ! CL 
            CL_newt = 2.0*sign_a*sa*sa*ca
            exp_pos = exp( lambda_CL*(alpha-alpha_0CL+alpha_sCL))
            exp_neg = exp(-lambda_CL*(alpha-alpha_0CL-alpha_sCL))
            sigma = (1.0 + exp_neg + exp_pos)/((1.0 + exp_neg)*(1.0 + exp_pos))
            CL = (1.0-sigma)*CL + sigma*CL_newt
            ! CD 
            CD_newt = 2.0*sin(abs(alpha))**3
            exp_pos = exp( lambda_CD*(alpha-alpha_0CD+alpha_sCD))
            exp_neg = exp(-lambda_CD*(alpha-alpha_0CD-alpha_sCD))
            sigma = (1.0 + exp_neg + exp_pos)/((1.0 + exp_neg)*(1.0 + exp_pos))
            CD = (1.0-sigma)*CD + sigma*CD_newt
            ! Cm
            Cm_newt = Cmmin*sign_a*sa*sa
            exp_pos = exp( lambda_Cm*(alpha-alpha_0Cm+alpha_sCm))
            exp_neg = exp(-lambda_Cm*(alpha-alpha_0Cm-alpha_sCm))
            sigma = (1.0 + exp_neg + exp_pos)/((1.0 + exp_neg)*(1.0 + exp_pos))
            Cm = (1.0-sigma)*Cm + sigma*Cm_newt
        end if 

        FM(1) = CL*sa-CS*ca*sb-CD*ca*cb
        FM(2) = CS*cb-CD*sb
        FM(3) = -CL*ca-CS*sa*sb-CD*sa*cb
        FM(4) = lat_ref*Cll
        FM(5) = long_ref*Cm
        FM(6) = lat_ref*Cn
        FM = 0.5*rho*V**2*sref*FM
        FM(1) = FM(1) + tau*Thrust0*(rho/rho0)**Ta
        FM(4:6) = FM(4:6) + cross_product_3D(aero_ref_location, FM(1:3)) ! body-fixed axis
        if (rk4_verbose) then
            write(*,*) "FM"
            write(*,*) FM
        end if
    end subroutine pseudo_aero

    subroutine print_aero_table()
        implicit none 
        integer :: i, iunit 
        real :: alpha, beta, states(13)
        real :: N_force, Y_force, A_force 
        real :: controls_original(4)
        real :: ca, cb, sa, sb 
        real :: CL, CD, Cm 
        real :: Z, T, P, rho, a, mu, const 
        write(*,*) "Starting to print aero tables"
        write(*,*) "Starting to print aero tables"
        write(*,*) "Starting to print aero tables"
        write(*,*) "Starting to print aero tables"

        call std_atm_English(Alt_initial, Z, T, P, rho, a, mu)
        const = 0.5*rho*V_initial**2*sref 

        controls_original = controls
        controls = 0.0
        states = 0.0

        open(newunit=iunit, file='aerotable.csv', status = 'REPLACE')
        write(iunit,*) 'alpha[deg],CL,CD,Cm'
        beta = 0.0
        do i=-180,180,1
            alpha = real(i)*PI/180.0
            beta = 0.0

            states(1) = V_initial*cos(alpha)
            states(2) = V_initial*sin(beta)
            states(3) = V_initial*sin(alpha)
            states(9) = -Alt_initial

            call pseudo_aero(states)
            A_force = -FM(1)
            Y_force =  FM(2)
            N_force = -FM(3)

            ca = cos(alpha)
            cb = cos(beta)
            sa = sin(alpha)
            sb = sin(beta)
            
            CL = N_force*ca - A_force*sa 
            CD = A_force*ca*cb - Y_force*sb + N_force*sa*cb 
            Cm = FM(5)
            
            CL = CL/const 
            CD = CD/const 
            Cm = Cm/const/long_ref
            write(iunit,*) alpha*180.0/PI, ',',CL,',',CD,',',Cm
        end do 
        close(iunit)
        controls = controls_original
    end subroutine print_aero_table


end module sim_m
