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
    end type pid_t 

    type controller_t
        type(pid_t), allocatable :: pid_controllers(:)
        integer :: num_pid 
        type(connection) :: pilot_conn 
        logical :: running = .false. 
    end type 

contains 

    subroutine controller_init(this, j_main) 
        implicit none 
        type(controller_t), intent(inout) :: this 
        type(json_value), pointer :: j_main, j_connections, j_pid, j_temp 
        integer :: i 
        logical :: found 

        write(*,*) 'Initializing Controller...'
        this%running = .true. 
        write(*,*) 'Initializing Connections...'
        call jsonx_get(j_main, 'connections', j_connections)
        call jsonx_get(j_connections, 'receive_pilot_commands', j_temp)
        call this%pilot_conn%init(j_temp, time=0.0)

        write(*,*) 'Initializing All PID Controllers...'
        call jsonx_get(j_main, 'PID', j_pid)
        this%num_pid = json_value_count(j_pid)
        allocate(this%pid_controllers(this%num_pid))

        do i =1, this%num_pid 
            call json_value_get(j_pid, i, j_temp)
            call pid_init(this%pid_controllers(i),j_temp)
        end do 

        write(*,*) 'Controller Initialization Sucessfully Completed'
    end subroutine controller_init


    function controller_update(this, states, time) result(ans)
        implicit none 
        type(controller_t), intent(inout) :: this 
        real, intent(in) :: states(21), time 
        real :: ans(4) ! 4 commands (aileron, elevator, rudder, throttle)
        real :: omega_command(3), omega_actual(3)
        real :: Z, Temp, P, rho, a, mu, dyp 

        ans = [0.0, -9.086019165449*PI/180.0,0.0,0.070182002357] ! turn off after testing of aileron is completed!!!!

        omega_actual = states(4:6)
        call std_atm_English(-states(9),Z,Temp,P,rho,a,mu)
        dyp = 0.5*rho*(states(1)**2 + states(2)**2 + states(3)**2) ! current dynamic pressure

        omega_command = this%pilot_conn%recv([time], time)
        ans(1) = pid_get_command(this%pid_controllers(1), omega_command(1),omega_actual(1),time,dyp)
        ans(2) = pid_get_command(this%pid_controllers(2), omega_command(2),omega_actual(2),time,dyp)
        ans(3) = pid_get_command(this%pid_controllers(3), omega_command(3),omega_actual(3),time,dyp)
    end function controller_update

    subroutine pid_init(this, j_pid)
        implicit none 
        type(pid_t), intent(inout) :: this 
        type(json_value), pointer :: j_pid 

        this%name = j_pid%name 
        write(*,*) ' Initializing PID for', this%name 
        call jsonx_get(j_pid, 'update_rate[hz]', this%update_rate)
        call jsonx_get(j_pid, 'KP', this%KP)
        call jsonx_get(j_pid, 'KI', this%KI)
        call jsonx_get(j_pid, 'KD', this%KD)

        this%prev_error = 0.0
        this%error_int = 0.0
        this%prev_ans = 0.0 
        this%prev_time = -1.0

        call jsonx_get(j_pid, 'units', this%units,'none') ! defaults to none 
        if(this%units == 'deg') then 
            this%display_units = 180.0/PI 
        else 
            this%display_units = 1.0 
        end if 
        call jsonx_get(j_pid, 'limits', this%limit,0.0,2)
        write(*,*) this%limit(:)
        this%limit(:) = this%limit(:)/this%display_units
        write(*,*) this%limit(:)
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
            this%prev_ans = ans 
            return 
        end if 

        if(time-this%prev_time >= 1.0/this%update_rate - TOLERANCE) then ! in this case, update the controller
            dt = time - this%prev_time 
            ! Error component 
            this%error = commanded - actual 
            ! Integrator with clamping 
            if ((this%prev_ans > this%limit(1)) .and. (this%prev_ans < this%limit(2))) then 
                this%error_int = this%error_int + 0.5*(this%prev_error + this%error)*dt
            else 
                write(*,*) this%name, ' PID controller saturated at ',&
                 this%prev_ans*this%display_units,'. Using Integrator clamping'
            end if 
            ! Derivative 
            if(dt>TOLERANCE) this%error_deriv = (this%error-this%prev_error)/dt 
            ! Final answer with proportional, integral, and derivative
            ans = (this%KP*this%error + this%KI*this%error_int + this%KD*this%error_deriv)/dyp 

            this%prev_error = this%error 
            this%prev_time = time 
            this%prev_ans = ans 
        else 
            ans = this%prev_ans 
        end if 
    end function pid_get_command


end module controller_m 