module adams_m
    implicit none
    real, parameter :: PI = 3.141592653589793
    real, parameter :: identity(3,3) = reshape([ 1.0, 0.0, 0.0, &
                                             0.0, 1.0, 0.0, &
                                             0.0, 0.0, 1.0 ], [3,3])
    real, parameter :: TOLERANCE = 1.0e-12
    real, parameter :: g_ssl = 9.80665 ! [m*s^-2] this is the gravitational acceleration at sea level
    real, parameter :: g_ssl_English = g_ssl/0.3048 ! [ft*s^-2] this is the gravitational acceleration at sea level
    real, parameter :: R_EZ = 6356766 ! [m] this is the earth radius for Geopotential altitude
    real, parameter :: R_EZ_English = R_EZ/0.3048 ! [ft]
    real, parameter :: R_E = 6366707.01949371  ! [m] this is the mean radius of earth at sea level
    real, parameter :: R_E_English = R_E/0.3048 ! [ft]
    real, parameter :: R = 287.0528  ! [J/(kg*K)] Gas constant for dry air
    real, parameter :: gamma = 1.4 ! Heat capacity ratio
    real, parameter :: sqrt_gammaR = sqrt(gamma*R) ! speeds up SI atmosphere test
    real, parameter :: mu_s = 1.716e-5 ! [kg/(m*s)] reference viscosity at T_s = 273.15 K
    real, parameter :: T_s = 273.15 ! [K] reference temperature for Sutherland's formula
    real, parameter :: K_s = 110.4 !  Sutherland's constant
    !!!Parameters for STD Atm in metric!!!
    ! real, parameter :: z0 = 0.0, z_11 = 11000.0, z_20 = 20000.0, &
    !                 z_32 = 32000.0, z_47 = 47000.0, z_52 = 52000.0, &
    !                 z_61 = 61000.0, z_79 = 79000.0
    real, parameter :: T0 = 288.150, T_11 = 216.650, T_20 = 216.650, &
                T_32 = 228.650, T_47 = 270.650, T_52 = 270.650, &
                T_61 = 252.650, T_79 = 180.650
    real, parameter :: T_prime0 = -6.5e-3, T_prime11 = 0.0, T_prime20 = 1.0e-3, &
                T_prime32 = 2.8e-3, T_prime47 = 0.0, T_prime52 = -2.0e-3, &
                T_prime61 = -4.0e-3, T_prime79 = 0.0
    real, parameter :: P0 = 1.01325e+05, P_11 = 2.26320318222212e+04, P_20 = 5.47487352827083e+03, &
                P_32 =  8.68014769086723e+02, P_47 =  1.10905588989225e+02, P_52 =  5.90005242789244e+01, &
                P_61 =  1.82099249050177e+01, P_79 = 1.03770045489203
    real, parameter :: invR_T0 = -g_ssl/(R*T_prime0), invR_T11 = -g_ssl/(R*T_11), invR_T20 = -g_ssl/(R*T_prime20), &
                       invR_T32 = -g_ssl/(R*T_prime32), invR_T47 = -g_ssl/(R*T_47), invR_T52 = -g_ssl/(R*T_prime52), & 
                       invR_T61 = -g_ssl/(R*T_prime61), invR_T79 = -g_ssl/(R*T_79)

contains

!!! ATTITUDE DESCRIPTORS (EULER, QUATERNIONS) !!!
function quat_mult(q1, q2) result(ans)
    implicit none
    real, intent(in) :: q1(4), q2(4)
    real :: ans(4)
    real :: e01, ex1, ey1, ez1
    real :: e02, ex2, ey2, ez2
    e01 = q1(1) 
    ex1 = q1(2) 
    ey1 = q1(3)
    ez1 = q1(4)
    e02 = q2(1) 
    ex2 = q2(2) 
    ey2 = q2(3)
    ez2 = q2(4)
    ans(1) = e01*e02 - ex1*ex2 - ey1*ey2 -ez1*ez2
    ans(2) = e01*ex2 + ex1*e02 + ey1*ez2 -ez1*ey2
    ans(3) = e01*ey2 - ex1*ez2 + ey1*e02 + ez1*ex2
    ans(4) = e01*ez2 + ex1*ey2 - ey1*ex2 + ez1*e02
end function quat_mult

function quat_base_to_dependent(vec, quat) result(v2)
    ! Eq. (1.5.4 in the book which uses the quaternion product twice to get Eq. 1.5.7)
    implicit none
    real, intent(in) :: vec(3), quat(4)
    real :: v2(3)
    real :: t0, tx, ty, tz
    real :: e0, ex, ey, ez
    real :: v0, vx, vy, vz
    v0 = 0
    vx = vec(1)
    vy = vec(2)
    vz = vec(3)
    e0 = quat(1)
    ex = quat(2)
    ey = quat(3)
    ez = quat(4)
    t0 = -vx*ex - vy* ey - vz*ez 
    tx =  vx*e0 + vy* ez - vz* ey 
    ty =  -vx* ez + vy* e0 + vz* ex 
    tz = vx* ey - vy* ex + vz* e0
    v2(1) = e0*tx - ex*t0 - ey*tz + ez*ty 
    v2(2) = e0*ty + ex*tz - ey*t0 - ez*tx 
    v2(3) = e0*tz - ex*ty + ey*tx - ez*t0
end function quat_base_to_dependent

function quat_dependent_to_base(v2, quat) result(v1)
    implicit none
    real, intent(in) :: v2(3), quat(4)
    real :: v1(3)
    real :: t0, tx, ty, tz
    real :: e0, ex, ey, ez
    real :: v0, vx, vy, vz
    v0 = 0
    vx = v2(1)
    vy = v2(2)
    vz = v2(3)
    e0 = quat(1)
    ex = quat(2)
    ey = quat(3)
    ez = quat(4)
    t0 = vx*ex + vy*ey + vz*ez 
    tx =  vx*e0 - vy*ez + vz*ey 
    ty =  vx*ez + vy*e0 - vz*ex 
    tz = -vx*ey + vy*ex + vz*e0
    v1(1) = e0*tx + ex*t0 + ey*tz - ez*ty 
    v1(2) = e0*ty - ex*tz + ey*t0 + ez*tx 
    v1(3) = e0*tz + ex*ty - ey*tx + ez*t0
end function quat_dependent_to_base

subroutine quat_norm(quat)
    ! normalizes a quaternion such that it becomes a unit quaternion. 
    ! Done the EXACT way a unit vector is the vector divided by its magnitude
    implicit none
    real, intent(inout) :: quat(4)
    real :: magnitude 
    magnitude = sqrt(quat(1)**2 + quat(2)**2 + quat(3)**2 + quat(4)**2)
    quat = quat/magnitude
end subroutine quat_norm

function euler_to_quat(euler) result(quat)
    ! Eq. 1.6.2 in the book. 
    real, intent(in) :: euler(3) 
    real :: quat(4)
    real :: phi, theta, psi 
    real :: cos_phi, cos_theta, cos_psi
    real :: sin_phi, sin_theta, sin_psi
    real :: cos_phi_cos_theta, sin_phi_sin_theta, sin_phi_cos_theta, cos_phi_sin_theta
    real :: hphi, htheta, hpsi
    phi = euler(1)
    theta = euler(2)
    psi = euler(3)
    hphi = 0.5*phi
    htheta = 0.5*theta
    hpsi = 0.5*psi
    cos_phi = cos(hphi)
    cos_theta = cos(htheta)
    cos_psi = cos(hpsi)
    sin_phi = sin(hphi)
    sin_theta = sin(htheta)
    sin_psi = sin(hpsi)
    cos_phi_cos_theta = cos_phi*cos_theta
    sin_phi_sin_theta = sin_phi*sin_theta
    sin_phi_cos_theta = sin_phi*cos_theta
    cos_phi_sin_theta = cos_phi*sin_theta
    quat(1) = cos_phi_cos_theta*cos_psi + sin_phi_sin_theta*sin_psi
    quat(2) = sin_phi_cos_theta*cos_psi - cos_phi_sin_theta*sin_psi
    quat(3) = cos_phi_sin_theta*cos_psi + sin_phi_cos_theta*sin_psi
    quat(4) = cos_phi_cos_theta*sin_psi - sin_phi_sin_theta*cos_psi
end function euler_to_quat

function quat_to_euler(quat) result(euler)
    ! Eq. 1.6.3 in the book. 
    real, intent(in) :: quat(4) 
    real :: euler(3)
    real :: cos_pi_4, arcsin_quat_cos_pi_4
    real :: e0, ex, ey, ez
    e0 = quat(1)
    ex = quat(2)
    ey = quat(3)
    ez = quat(4)
    cos_pi_4 = ex/cos(PI*0.25)
    arcsin_quat_cos_pi_4 = asin(cos_pi_4)
    if (e0*ey-ex*ez == 0.5) then
        euler(1) = 2*arcsin_quat_cos_pi_4
        euler(2) = 0.5*PI
        euler(3) = 0.0
    else if (e0*ey-ex*ez == 0.5) then
        euler(1) = 2*arcsin_quat_cos_pi_4
        euler(2) = -0.5*PI
        euler(3) = 0.0
    else 
        euler(1) = atan2(2*(e0*ex + ey*ez), (e0**2 + ez**2 -ex**2 -ey**2))
        euler(2) = asin(2*(e0*ey-ex*ez))
        euler(3) = atan2(2*(e0*ez + ex*ey), (e0**2 + ex**2 -ey**2 -ez**2))
        if (euler(3) < 0) then 
            euler(3) = euler(3) + 2*PI 
        end if 
    end if
end function quat_to_euler

!!! STANDARD ATMOSPHERE !!!
function gravity_SI(H) result(g) ! 3.13.1
    ! H is the Geometric height above sea-level in meters. 
    real, intent(in) :: H
    real :: g
    g = g_ssl*(R_EZ/(R_EZ+H))**2
end function gravity_SI

function gravity_English(H) result(g) ! 3.13.2
    ! H is the Geometric height above sea-level in ft. 
    real, intent(in) :: H
    real :: g
    g = g_ssl_English*((R_EZ_English)/(R_EZ_English+H))**2
end function gravity_English

subroutine std_atm_SI(H, Z, T, P, rho, a, mu)
    ! Z = Geopotential Alt (Eq 3.2.2), T = temp (from Eq 3.2.3 and table 3.2.1), P = pressure (Eq. 3.2.7 and table 3.2.1), rho = density (Eq. 3.2.8), a = sonic velocity (Eq, 3.2.9)
    real, intent(in) :: H
    real, intent(inout) ::  Z, T, P, rho, a, mu
    ! first calculate Geopotential Alt
    Z = R_EZ*H/(R_EZ+H)
    Z = max(Z, 0.0)
    ! now calculate the temperature and pressure at (Z) depending on what range H is in 
    if (Z < 11000.0) then
        T = T0 + T_prime0*Z
        P = P0*((T)/(T0))**(invR_T0)
    else if (Z >= 11000.0 .and. Z < 20000.0) then 
        T = T_11
        P = P_11*exp(invR_T11*(Z-11000.0))
    else if (Z >= 20000.0 .and. Z < 32000.0) then
        T = T_20 + T_prime20*(Z-20000.0)
        P = P_20*((T)/(T_20))**(invR_T20)
    else if (Z >= 32000.0 .and. Z < 47000.0) then
        T = T_32 + T_prime32*(Z-32000.0)
        P = P_32*((T)/(T_32))**(invR_T32)
    else if (Z >= 47000.0 .and. Z < 52000.0) then 
        T = T_47
        P = P_47*exp(invR_T47*(Z-47000.0))
    else if (Z >= 52000.0 .and. Z < 61000.0) then
        T = T_52 + T_prime52*(Z-52000.0)
        P = P_52*((T)/(T_52))**(invR_T52)
    else if (Z >= 61000.0 .and. Z < 79000.0) then
        T = T_61 + T_prime61*(Z-61000.0)
        P = P_61*((T)/(T_61))**(invR_T61)
    else  
        T = T_79 
        P = P_79*exp(invR_T79*(Z-79000.0))
    end if
    ! calculate density from pressure, R (gas const for dry air), and Temperature
    rho = P/(R*T)
    ! calculate speed of sound from gamma (heat capacity ratio), R (gas const dry air), and Temp
    a = sqrt_gammaR*sqrt(T)
    ! calculate dynamic viscosity using Sutherland's formula (Eq. 3.2.10)
    mu = mu_s*(T_s+K_s)/(T+K_s)*(T/T_s)**1.5
end subroutine std_atm_SI

subroutine std_atm_English(H, Z, T, P, rho, a, mu)
    real, intent(in) :: H
    real, intent(inout) ::  Z, T, P, rho, a, mu
    real :: si_H
    !! converts H into SI units, calls std_atm_SI, converts H, Z, T< P, rho, and a back to english units
    si_H = H*0.3048 ! converting ft to meters for the wrapper
    call std_atm_SI(si_H, Z, T, P, rho, a, mu)
    !! convert Z, T, P, rho, and a back into English units
    Z = Z/0.3048
    T = T*1.8
    P = P/47.880258
    rho = rho/515.379
    a = a/0.3048
    mu = mu/47.880258
end subroutine std_atm_English

!!!! TURBULENCE STUFF BELOW 
!---------------------------------------------------------------------------------
function rand_normal() result(z)
    !! Return a single normally distributed random number (mean=0, std=1)
    real :: z
    real, save :: z2 = 0.0
    logical, save :: has_saved = .false.
    real :: u1, u2, r, theta
    
    if (has_saved) then
        !! Return the saved value from last time
        z = z2
        has_saved = .false.
    else
        !! Generate two uniforms
        call random_number(u1)
        call random_number(u2)
        
        !! Guard against u1 = 0
        if (u1 <= 0.0) u1 = 1.0e-12
        if (u1 >  1.0) u1 = 1.0
        
        !! Box-Muller transform
        r = sqrt(-2.0*log(u1))
        theta = 2.0*pi*u2
        z  = r * cos(theta)
        z2 = r * sin(theta)
        
        has_saved = .true.
    end if
end function rand_normal

!---------------------------------------------------------------------------------
function   interpolate_1D(xa, ya, x, saturate) result(y)
    real, intent(in) :: xa(:), ya(:), x
    logical, intent(in), optional :: saturate
    real :: y
    integer :: i, N
    logical :: sat
    N = size(xa)
    if (N /= size(ya)) then
        write(*,*) 'Error interpolating 1D array: sizes do not match!'
        stop
    end if
    if (present(saturate)) then
        sat = saturate
    else
        sat = .true.
    end if
    if (sat) then
        if (x <= xa(1)) then
            y = ya(1)
            return
        end if
        do i=2, N
            if (xa(i-1) <= x .and. x <= xa(i)) then
                y = interpolate(xa(i-1), xa(i), x, ya(i-1), ya(i))
                return
            end if
        end do
        y = ya(N)
    else
        do i=2, N
            if (xa(i-1) <= x .and. x <= xa(i)) then
                y = interpolate(xa(i-1), xa(i), x, ya(i-1), ya(i))
                return
            end if
        end do
        write(*,*) 'Error interpolating 1D array: value outside of the range'
        stop
    end if
end function interpolate_1D

!---------------------------------------------------------------------------------
function interpolate(x1, x2, x, y1, y2) result(y)
    real, intent(in) :: x1, x2, x, y1, y2
    real :: y
    y = (y2-y1) / (x2-x1) * (x-x1) + y1
end function interpolate

!---------------------------------------------------------------------------------
subroutine psd(x, dt, psd_norm, filename)
    implicit none
    real, intent(in) :: x(:), dt
    real, intent(out), optional :: psd_norm(:,:)
    character(len=*), intent(in), optional :: filename
    integer :: n, k, i, iunit
    real :: fs, f, mean, stdev
    real, allocatable :: Pxx(:), loc_psd_norm(:,:)
    complex :: eye, temp

    n = size(x)
    if (mod(n,2) /= 0) error stop "periodogram: N must be even for N/2 bin."
    allocate(Pxx(n/2+1))
    allocate(loc_psd_norm(n/2+1,2))
    eye = cmplx(0.0, 1.0)

    ! compute mean and stdev
    mean = 0.0
    stdev = 0.0
    do i=1,n
        mean = mean + x(i)/n
    end do
    do i=1,n
        stdev = stdev + (x(i)-mean)**2
    end do
    stdev = sqrt(stdev/(n-1))

    if(present(filename)) then
        open(newunit=iunit, file=filename, status="replace", action="write")
        write(iunit,'(A)') 'fs[Hz],PSD[ft^2/Hz],Normalized_PSD,Mean,STDEV'
        write(*,*) trim(filename),'   Mean = ',mean,'   Std. Dev = ',stdev
    end if

    Pxx(:) = 0.0
    fs = 1.0/dt
    do k=0,n/2
        f = real(k)*fs/real(n)
        temp = cmplx(0.0, 0.0)
        do i=1,n
            temp = temp + (x(i)-mean)*exp(-eye*2.0*PI*real(i-1)*real(k)/real(n))
        end do
        Pxx(k+1) = 1/fs/real(n)*abs(temp)**2
        if (k /= 0 .and. k /= n/2) Pxx(k+1) = 2.0*Pxx(k+1)
        loc_psd_norm(k+1,1) = f
        loc_psd_norm(k+1,2) = Pxx(k+1)/stdev**2
        if(present(filename)) then
            if(k==0) then
                write(iunit,'(ES20.12,",",ES20.12,",",ES20.12,",",ES20.12,",",ES20.12)') f, Pxx(k+1), loc_psd_norm(k+1,2),mean,stdev
            else
                write(iunit,'(ES20.12,",",ES20.12,",",ES20.12)') f, Pxx(k+1), loc_psd_norm(k+1,2)
            end if
        end if
    end do

    if(present(filename)) close(iunit)
    if(present(psd_norm)) psd_norm = loc_psd_norm
end subroutine psd

! subroutine psd(x, dt, psd_norm, filename)
!     implicit none
!     real, intent(in) :: x(:), dt
!     real, intent(out), optional :: psd_norm(:,:)
!     character(len=*), intent(in), optional :: filename
!     integer :: n, k, i, iunit
!     real :: fs, f, mean, stdev
!     real, allocatable :: Pxx(:), loc_psd_norm(:,:)
!     complex, allocatable :: cx(:)

!     n = size(x)
!     if (mod(n,2) /= 0) error stop "periodogram: N must be even."

!     allocate(cx(n))                         ! <-- NO PADDING
!     allocate(Pxx(n/2+1))
!     allocate(loc_psd_norm(n/2+1,2))

!     ! Mean
!     mean = 0.0
!     do i=1,n
!         mean = mean + x(i)/n
!     end do

!     ! Std dev
!     stdev = 0.0
!     do i=1,n
!         stdev = stdev + (x(i)-mean)**2
!     end do
!     stdev = sqrt(stdev/(n-1))

!     if(present(filename)) then
!         open(newunit=iunit, file=filename, status="replace", action="write")
!         write(iunit,'(A)') 'fs[Hz],PSD[ft^2/Hz],Normalized_PSD,Mean,STDEV'
!         write(*,*) trim(filename),'   Mean = ',mean,'   Std. Dev = ',stdev
!     end if

!     ! Build input (IDENTICAL to DFT definition)
!     do i = 1, n
!         cx(i) = cmplx(x(i) - mean, 0.0)
!     end do

!     call fft(cx)   ! length = n

!     fs = 1.0/dt

!     do k = 0, n/2
!         f = real(k)*fs/real(n)

!         ! EXACT SAME normalization as your DFT
!         Pxx(k+1) = (1.0/fs/real(n)) * abs(cx(k+1))**2

!         if (k /= 0 .and. k /= n/2) Pxx(k+1) = 2.0*Pxx(k+1)

!         loc_psd_norm(k+1,1) = f
!         loc_psd_norm(k+1,2) = Pxx(k+1)/stdev**2

!         if(present(filename)) then
!             if(k==0) then
!                 write(iunit,'(ES20.12,",",ES20.12,",",ES20.12,",",ES20.12,",",ES20.12)') f, Pxx(k+1), loc_psd_norm(k+1,2),mean,stdev
!             else
!                 write(iunit,'(ES20.12,",",ES20.12,",",ES20.12)') f, Pxx(k+1), loc_psd_norm(k+1,2)
!             end if
!         end if
!     end do

!     if(present(filename)) close(iunit)
!     if(present(psd_norm)) psd_norm = loc_psd_norm
! end subroutine psd

! function next_power_of_2(n) result(m)
!         implicit none
!         integer, intent(in) :: n
!         integer :: m
!         m = 1
!         do while (m < n)
!             m = m * 2
!         end do
! end function next_power_of_2

! subroutine fft(x)
!     complex, intent(inout) :: x(:)
!     integer :: n, m, j, k, half, i
!     complex :: t, w, wn

!     n = size(x)

!     ! Bit-reversal permutation
!     j = 1
!     do i = 1, n-1
!         if (i < j) then
!             t = x(j); x(j) = x(i); x(i) = t
!         end if
!         m = n/2
!         do while (m >= 1 .and. j > m)
!             j = j - m; m = m/2
!         end do
!         j = j + m
!     end do

!     ! Cooley-Tukey butterfly stages
!     half = 1
!     do while (half < n)
!         wn = exp(cmplx(0.0, -2.0*PI/real(2*half)))
!         w  = cmplx(1.0, 0.0)
!         do k = 1, half
!             do j = k, n, 2*half
!                 t          = w * x(j + half)
!                 x(j + half) = x(j) - t
!                 x(j)       = x(j) + t
!             end do
!             w = w * wn
!         end do
!         half = half * 2
!     end do
! end subroutine fft

!---------------------------------------------------------------------------------
subroutine test_rand_normal()
    implicit none
    integer, parameter :: n = 100000
    integer, parameter :: nbins = 100
    real, allocatable :: x(:), hist(:)
    real :: mean, sigma, xmin, xmax, dx, binval
    integer :: i, b, iunit
    character(len=:), allocatable :: filename  

    write(*,*) '   Testing function rand_normal():'

    allocate(x(n))
    allocate(hist(nbins))
    do i = 1, n
        ! call random_number(x(i)) ! this is a uniform distribution from 0 to 1 with a mean of 0.5
        x(i) = rand_normal() ! this is a normal distribution with a mean of 0 and a std dev of 1
    end do

    mean  = sum(x) / n
    sigma = sqrt( sum((x-mean)**2) / (n-1) )

    write(*,*) '     Mean  = ', mean
    write(*,*) '     Std   = ', sigma
    write(*,*) '     Expected mean error ~ ', 1.0/sqrt(real(n))
    write(*,*) '     Expected std error  ~ ', 1.0/sqrt(2.0*n)

    ! Build histogram
    xmin = -5.0; xmax = 5.0
    dx = (xmax-xmin)/nbins
    hist = 0.0

    do i = 1, n
        if (x(i) >= xmin .and. x(i) < xmax) then
            b = int((x(i)-xmin)/dx) + 1
            hist(b) = hist(b) + 1
        end if
    end do

    hist = hist / (n*dx)   ! normalize to PDF

    filename = 'rand_normal_test.csv'
    open(newunit=iunit, file=filename, status="replace", action="write")
    write(iunit,'(A)') 'bin_number,bin_midpoint_value,normalized_histogram,analytic_solution'
    do i=1,nbins
        binval = xmin+dx*i-0.5*dx
        write(iunit,*) i,',',binval,',',hist(i),',',exp(-0.5*(binval**2))/sqrt(2.0*PI)
    end do
    write(*,*) '     Histogram written to file: ',trim(filename)
    write(*,*)
    close(iunit)

    ! Compute PSD
    ! call psd(x(:),0.01,filename='rand_normal_psd.csv')

    deallocate(x)
    deallocate(hist)
end subroutine test_rand_normal

end module adams_m
