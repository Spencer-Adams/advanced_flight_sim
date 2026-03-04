module controller_m 
    use adams_m
    use jsonx_m
    use linalg_mod
    use micro_time_m
    use connection_m
    implicit none 

    type pid_t
        character(len=:), allocatable :: name 
        real :: KP, KI, KD 
        real :: error, prev_error, error_int, error_deriv 
        real :: prev_time, prev_ans, update_rate 
        real, allocatable :: limit(:)
        character(len=:), allocatable :: units 
        real :: display_units = 1.0 
        logical :: dyp_schedule 
    end type pid_t 

    type controller_t
        type(pid_t) :: p_da, q_de, r_dr, bank_p, gamma_q, V_tau
        type(connection) :: pilot_conn 
        logical :: running = .false. 
    end type 

contains 

    subroutine controller_init(this, j_main) 
        implicit none 
        type(controller_t), intent(inout) :: this 
        type(json_value), pointer :: j_main, j_connections, j_pid, j_temp 
        logical :: found 

        write(*,*) 'Initializing Controller...'
        this%running = .true. 
        write(*,*) 'Initializing Connections...'
        call jsonx_get(j_main, 'connections', j_connections)
        call jsonx_get(j_connections, 'receive_pilot_commands', j_temp)
        call this%pilot_conn%init(j_temp, time=0.0)

        write(*,*) 'Initializing All PID Controllers...'
        call jsonx_get(j_main, 'PID', j_pid)

        ! p to aileron 
        call jsonx_get(j_pid, 'p->aileron', j_temp)
        call pid_init(this%p_da, j_temp)

        ! q to elevator 
        call jsonx_get(j_pid, 'q->elevator', j_temp)
        call pid_init(this%q_de, j_temp)

        ! r to rudder 
        call jsonx_get(j_pid, 'r->rudder', j_temp)
        call pid_init(this%r_dr, j_temp)

        ! bank to p 
        call jsonx_get(j_pid, 'bank->p', j_temp)
        call pid_init(this%bank_p, j_temp)

        ! climb to q 
        call jsonx_get(j_pid, 'gamma->q', j_temp)
        call pid_init(this%gamma_q, j_temp)

        ! V to throttle 
        call jsonx_get(j_pid, 'V->throttle', j_temp)
        call pid_init(this%V_tau, j_temp)

        write(*,*) "Controller Initialization Complete."
    end subroutine controller_init

    function controller_update(this, states, time) result(ans)
        implicit none 
        type(controller_t), intent(inout) :: this 
        real, intent(in) :: states(21), time 
        real :: ans(4) ! 4 commands (aileron, elevator, rudder, throttle)
        real :: pilot_command(3), g, gamma 
        real :: bank_sp, gamma_sp, V_sp
        real :: u, v, w, p, q, r, eul(3), sp, cp, st, ct, Vmag
        real :: p_sp, q_sp, r_sp 
        real :: Z_temp, T_temp, P_temp, a_temp, mu_temp 
        real :: rho, dyp 

        u = states(1)
        v = states(2)
        w = states(3)
        p = states(4)
        q = states(5)
        r = states(6)

        eul = quat_to_euler(states(10:13))
        sp = sin(eul(1)) 
        st = sin(eul(2)) 
        cp = cos(eul(1)) 
        ct = cos(eul(2))  
        Vmag = sqrt(u**2 + v**2 + w**2)
        gamma = asin((u*st-(v*sp + w*cp)*ct)/Vmag)

        call std_atm_English(-states(9), Z_temp, T_temp, P_temp, rho, a_temp, mu_temp)
        dyp = 0.5*rho*Vmag**2 
        g = gravity_English(-states(9))

        pilot_command = this%pilot_conn%recv([time],time)
        bank_sp = pilot_command(1)*PI/180.0
        gamma_sp = pilot_command(2)*PI/180.0
        V_sp = pilot_command(3)

        p_sp = pid_get_command(this%bank_p, bank_sp, eul(1), time, dyp)
        q_sp = pid_get_command(this%gamma_q, gamma_sp, gamma, time, dyp)
        r_sp = (g*sp*ct + p*w)/u ! no gravity relief 

        !! Do any pilot command 
        ans(1) = pid_get_command(this%p_da, p_sp, p, time, dyp)
        ans(2) = pid_get_command(this%q_de, q_sp, q, time, dyp)
        ans(3) = pid_get_command(this%r_dr, r_sp, r, time, dyp)
        ans(4) = pid_get_command(this%V_tau, V_sp, Vmag, time, dyp)

        !! 12.6.1 
        ! pilot_12.6.1.csv 
        ! ans(1) = pid_get_command(this%p_da, bank_sp, p, time, dyp) ! aileron
        ! write(*,*) "p_sp = ", bank_sp, ", p = ", p 
        ! ans(2) = -9.086019165449*PI/180.0 ! elevator 
        ! ans(3) = 0.0 ! rudder
        ! ans(4) = 0.070182002357 ! throttle 
        
        !! 12.6.2 
        ! ans(1) = pid_get_command(this%p_da,0.0,p,time,dyp)
        ! ans(2) = pid_get_command(this%q_de,0.0,q,time,dyp)
        ! ans(3) = pid_get_command(this%r_dr,0.0,r,time,dyp)
        ! ans(4) = 0.070182002357 ! throttle 
    end function controller_update

    subroutine pid_init(this, j_pid)
        implicit none 
        type(pid_t), intent(inout) :: this 
        type(json_value), pointer :: j_pid 

        this%name = j_pid%name 
        write(*,*) "Name = ", this%name 
        call jsonx_get(j_pid, 'update_rate[hz]', this%update_rate)
        call jsonx_get(j_pid, 'kp', this%KP)
        call jsonx_get(j_pid, 'ki', this%KI)
        call jsonx_get(j_pid, 'kd', this%KD)

        this%prev_error = 0.0
        this%error_int = 0.0
        this%prev_ans = 0.0 
        this%prev_time = -1.0

        call jsonx_get(j_pid, 'units', this%units,'none') ! defaults to none 
        call jsonx_get(j_pid, 'dynamic_pressure_schedule', this%dyp_schedule) ! defaults to false 
        if(this%units == 'deg') then 
            this%display_units = 180.0/PI 
        else 
            this%display_units = 1.0 
        end if 
        call jsonx_get(j_pid, 'limits', this%limit,0.0,2)
        write(*,*) "limits not divided by display units: ", this%limit(:)
        this%limit(:) = this%limit(:)/this%display_units
        write(*,*) "limits divided by display units: ", this%limit(:)
    end subroutine pid_init

    function pid_get_command(this, commanded,actual,time,dyp) result(ans)
        implicit none 
        type(pid_t), intent(inout) :: this 
        real, intent(in) :: commanded, actual, time, dyp 
        real :: dt, ans 

        if (this%prev_time < 0.0) then 
            this%prev_time = time 
            this%prev_error = commanded - actual 
            this%error_int = 0.0
            this%error_deriv = 0.0
            ans = this%KP*this%prev_error
            if(this%dyp_schedule) ans = ans/dyp 
            this%prev_ans = ans 
            return 
        end if 

        if(time-this%prev_time >= 1.0/this%update_rate - TOLERANCE) then ! update the controller
            dt = time - this%prev_time 
            ! Error component 
            this%error = commanded - actual 
            ! if this%name == ""
            ! Integrator with clamping 
            if ((this%prev_ans > this%limit(1)) .and. (this%prev_ans < this%limit(2))) then 
                this%error_int = this%error_int + 0.5*(this%prev_error + this%error)*dt
            end if 
            ! Derivative 
            if(dt>TOLERANCE) this%error_deriv = (this%error-this%prev_error)/dt 
            ! Final answer with proportional, integral, and derivative
            ans = (this%KP*this%error + this%KI*this%error_int + this%KD*this%error_deriv)
            if(this%dyp_schedule) ans = ans/dyp 

            this%prev_error = this%error 
            this%prev_time = time 
            this%prev_ans = ans 
        else 
            ans = this%prev_ans 
        end if 
    end function pid_get_command


end module controller_m 