module vehicle_m
    use adams_m
    use jsonx_m
    use linalg_mod
    use micro_time_m
    use connection_m

    type vehicle_t
        type(json_value), pointer :: j_vehicle
        
        character(len=:), allocatable :: name
        character(len=:), allocatable :: type
        character(100) :: states_filename, rk4_filename

        logical :: run_physics
        logical :: save_states, rk4_verbose
        integer :: iunit_states, iunit_rk4, iunit_trim
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
        character(len=:), allocatable :: init_type 
        real, allocatable :: thrust_orientation(:)
        real :: Z_temp,T_temp,P_temp,a_temp,mu_temp

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
                write(this%iunit_states,*) "        t[s]                   u[ft/s]               v[ft/s]"// &
                        "                 w[ft/s]                p[rad/s]             q[rad/s]"  // &
                        "            r[rad/s]                x[ft]                  y[ft]" // &                   
                        "                  z[ft]                    e0                   ex" // &
                        "                      ey                   ez"
                write(*,*) '   - saving states to ', this%states_filename
            end if 

            write(*,*) '   - Initializing mass and inertia properties'
            call mass_inertia(this)
            write(*,*) '   - Initializing aerodynamics'
            call jsonx_get(this%j_vehicle, "aerodynamics.reference.area[ft^2]", this%sref)
            call jsonx_get(this%j_vehicle, "aerodynamics.reference.longitudinal_length[ft]", this%long_ref)
            call jsonx_get(this%j_vehicle, "aerodynamics.reference.lateral_length[ft]", this%lat_ref)
            call jsonx_get(this%j_vehicle, "aerodynamics.reference.relative_location[ft]", this%aero_ref_location,0.0,3)
            
            call jsonx_get(this%j_vehicle, "thrust.Thrust0[lbf]", this%Thrust0)
            call jsonx_get(this%j_vehicle, "thrust.Ta", this%Ta)
            call jsonx_get(this%j_vehicle, "thrust.location[ft]", this%thrust_location,0.0,3)
            call jsonx_get(this%j_vehicle, "thrust.orientation[deg]", thrust_orientation, 0.0,3)
            thrust_orientation = thrust_orientation*PI/180.0
            this%thrust_quat = euler_to_quat(thrust_orientation)
            if (this%type == 'arrow') then
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CL.alpha", this%CLa)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.L0", this%CDL0)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.CD.CL1_CL1", this%CDL2)
                call jsonx_get(this%j_vehicle, "aerodynamics.coefficients.Cl.0", this%Cll0)
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
            call jsonx_get(this%j_vehicle, "initial.airspeed[ft/sec]", this%init_V)
            call jsonx_get(this%j_vehicle, "initial.altitude[ft]", this%init_state(9))
            this%init_state(9) = - this%init_state(9)
            this%init_alt = this%init_state(9)
            call jsonx_get(this%j_vehicle, "initial.Euler_angles[deg]", this%init_eul,0.0,3)
            this%init_eul = this%init_eul*PI/180.0
            call jsonx_get(this%j_vehicle, "initial.type", init_type)

            if (init_type=="state") then
                call init_to_state(this)
            end if 
            this%init_state(10:13) = euler_to_quat(this%init_eul)
            this%state = this%init_state

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
        type(json_value), pointer :: j_initial, j_state 
        real :: alpha, beta
        write(*,*) ' Setting Trim. This is not ready in the code, so DO NOT trust results if it does work at all'
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
        call jsonx_get(this%j_vehicle, "mass.weight[lbf]", weight)
        call jsonx_get(this%j_vehicle, "mass.Ixx[slug-ft^2]", this%I(1,1))
        call jsonx_get(this%j_vehicle, "mass.Iyy[slug-ft^2]", this%I(2,2))
        call jsonx_get(this%j_vehicle, "mass.Izz[slug-ft^2]", this%I(3,3))
        call jsonx_get(this%j_vehicle, "mass.Ixy[slug-ft^2]", this%I(1,2),0.0)
        call jsonx_get(this%j_vehicle, "mass.Ixz[slug-ft^2]", this%I(3,1),0.0)
        call jsonx_get(this%j_vehicle, "mass.Iyz[slug-ft^2]", this%I(3,2),0.0)
        call jsonx_get(this%j_vehicle, "mass.hx[slug-ft^2/s]", this%hx) ! in the 3 by 3 h matrix in Eq. 5.4.6
        call jsonx_get(this%j_vehicle, "mass.hy[slug-ft^2/s]", this%hy)
        call jsonx_get(this%j_vehicle, "mass.hz[slug-ft^2/s]", this%hz)
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
        gravity = gravity_English(0.0)
        this%mass = weight/gravity
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
        real :: exp_pos, exp_neg, sigma
        real :: CL_newt, CD_newt, Cm_newt
        real :: CL1_blended, CD1_blended, Cm1_blended
        ahat = 0.0
        !!!! receive controls from python script here !!!!
        da = this%controls(1)
        de = this%controls(2)
        dr = this%controls(3)
        tau = this%controls(4)
        if (tau < 0.0) then 
            tau = 0.0
        else if (tau > 1.0) then 
            tau = 1.0
        end if 

        call std_atm_English(-y(9), Z, T, P, rho, a, mu)

        V = sqrt(y(1)**2 + y(2)**2 + y(3)**2)
        alpha = atan2(y(3), y(1)) ! Eq. 3.4.4
        beta = asin(y(2)/V) ! Eq. 3.4.5
        pbar = 0.5*y(4)*this%lat_ref/(V)
        qbar = 0.5*y(5)*this%long_ref/(V)
        rbar = 0.5*y(6)*this%lat_ref/(V)

        sa = sin(alpha)
        ca = cos(alpha)
        sb = sin(beta)
        cb = cos(beta)
        sign_a = sign(1.0, alpha)

        CL1 = this%CL0 +this%CLa*alpha
        CL = CL1 + this%CLqbar*qbar+this%CLahat*ahat + this%CLde*de
        CS = this%CSb*beta + (this%CSpbar+this%CSapbar*alpha)*pbar + this%CSrbar*rbar + this%CSda*da + this%CSdr*dr
        CD = this%CDL0 + this%CDL1*CL1 + this%CDL2*CL1**2 + this%CDS2*CS**2 + &
            (this%CDqbar + this%CDaqbar*alpha)*qbar + (this%CDde + this%CDade*alpha)*de + this%CDde2*de**2
        Cll = this%Clb*beta + this%Clpbar*pbar + (this%Clrbar + this%Clarbar*alpha)*rbar + this%Clda*da + this%Cldr*dr
        Cm = this%Cm0 + this%Cma*alpha + this%Cmqbar*qbar + this%Cmahat*ahat + this%Cmde*de 
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
        if (this%rk4_verbose) then
            write(*,*) "FM"
            write(*,*) this%FM
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
        call pseudo_aero(this, state)
        rot_and_inertia_temp = & 
        [this%FM(4) + dot_product(this%h_gyro(1,:),pqr_temp) + ((this%Iyyb-this%Izzb)*q*r-this%Iyzb*(q**2-r**2)-&
        this%Ixzb*p*q+this%Ixyb*p*r)-this%hdot_gyro(1),&
        this%FM(5) + dot_product(this%h_gyro(2,:),pqr_temp) + ((this%Izzb-this%Ixxb)*p*r-this%Ixzb*(r**2-p**2)-&
        this%Ixyb*q*r+this%Iyzb*p*q)-this%hdot_gyro(2),& 
        this%FM(6) + dot_product(this%h_gyro(3,:),pqr_temp) + ((this%Ixxb-this%Iyyb)*p*q-this%Ixyb*(p**2-q**2)-&
        this%Iyzb*p*r+this%Ixzb*q*r)-this%hdot_gyro(3)]
        res(1) = 1/this%mass * this%FM(1) + gravity * 2 *(ex*ez-ey*e0) + r*v - q*w ! udot body-fixed
        res(2) = 1/this%mass * this%FM(2) + gravity * 2 *(ey*ez+ex*e0) + p*w - r*u ! vdot body-fixed
        res(3) = 1/this%mass * this%FM(3) + gravity * (ez**2+e0**2-ex**2-ey**2) + q*u - p*v ! wdot body-fixed
        res(4:6) = matmul(this%Iinv, rot_and_inertia_temp) ! pdot, qdot, rdot body-fixed
        res(7:9) = quat_dependent_to_base((/u,v,w/), (/e0, ex, ey, ez/)) !!! + wind ! xdot ydot zdot earth-fixed
        res(10) = 0.5 * dot_product((/-ex, -ey, -ez/),pqr_temp) !e0
        res(11) = 0.5 * dot_product((/e0, -ez, ey/),pqr_temp) !ex
        res(12) = 0.5 * dot_product((/ez, e0, -ex/),pqr_temp) !ey
        res(13) = 0.5 * dot_product((/-ey, ex, e0/),pqr_temp) !ez
        if (this%rk4_verbose) then
            write(*,*) "state"
            write(*,*) state
            write(*,*) ""
            write(*,*) "res"
            write(*,'(14E22.13)') res
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

        y = this%state 
        y1 = runge_kutta(this,t,y,dt)
        call quat_norm(y1(10:13))
        this%state = y1 
        if (sqrt(y1(4)**2+y1(5)**2+y1(6)**2)/2/PI*dt>0.1) then 
            write(*,*) 'Warning, rotation rates large for Rk4. See Eq. 5.7.3 in the book'
        end if 
         if (this%save_states) then 
            call vehicle_write_state(this, t+dt,y1)
        end if
        ! if(rk4_verbose) then 
        ! end if
    end subroutine vehicle_tick_state

    subroutine vehicle_write_state(this, time, state)
        implicit none 
        type(vehicle_t) :: this
        real, intent(in) :: time, state(13)
        logical :: is_open
        inquire(file=this%states_filename, opened = is_open)
        if (is_open) then 
            write(*,*) 'output.txt is already open. THIS PROGRAM WILL NOT RUN IF THAT IS NOT CLOSED'
        else 
            write(this%iunit_states,'(14E22.13)') time,state(:)
        end if 
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
    end function runge_kutta

end module vehicle_m