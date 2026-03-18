module atmosphere_m 
    use adams_m 
    use jsonx_m 
    implicit none 

    type atmosphere_t 
        real, allocatable :: wind(:) 
        character(len=:), allocatable :: turb_model, turb_intensity 
        real :: wingspan, hstab_dist, vstab_dist 
        logical :: turb_repeatable 

        real :: light_hag(3) = (/2000.,8000.,17000./) ! hag means height above ground 
        real :: light_sig(3) = (/5.,5.,3./)
        real :: moderate_hag(3) = (/2000.,11000.,45000./)! hag means height above ground 
        real :: moderate_sig(3) = (/10.,10.,3./)
        real :: severe_hag(4) = (/2000.,4000.,20000.,80000./) ! hag means height above ground 
        real :: severe_sig(4) = (/15.,21.,21.,3./)

        real, allocatable :: turb_hag(:), turb_sig(:)
        real :: w_disturbs(20), v_disturbs(20)
        real :: etau_array(20), etav_array(20), etaw_array(20), etap_array(20) 
        real :: prev_turb(4), prev_xyz(3), prev_f, prev_g
        real :: Lu, Lv, Lw, Lb 
        real :: xff 
    end type atmosphere_t
contains 

    subroutine atmosphere_init(this, j_atmosphere)
        implicit none 
        type(atmosphere_t) :: this 
        type(json_value), pointer :: j_atmosphere, j_turb, j_sample 
        logical :: found
        integer :: i, n, seed 
        integer, allocatable :: seed_array(:) 

        write(*,*) 'Initializing Atmospheric Model...'
        call jsonx_get(j_atmosphere, 'constant_wind[ft/s]', this%wind, 0.0, 3)
        write(*,*) "     Constant Wind [ft/s] = ", this%wind(:)
        call json_get(j_atmosphere, 'turbulence', j_turb, found)
        if(found) then 
            call jsonx_get(j_turb, 'model', this%turb_model, 'none')
            if(this%turb_model .ne. 'none') then 
                call jsonx_get(j_turb, 'wingspan[ft]', this%wingspan)
                call jsonx_get(j_turb, 'hstab_distance[ft]', this%hstab_dist)
                call jsonx_get(j_turb, 'vstab_distance[ft]', this%vstab_dist)
                call jsonx_get(j_turb, 'intensity', this%turb_intensity)
                call jsonx_get(j_turb, 'repeatable', this%turb_repeatable)
                write(*,*) "       Turbulence Intensity = ", this%turb_intensity
                write(*,*) "       Turbulence Model = ", this%turb_model 
                ! random number generator 
                if(this%turb_repeatable) then 
                    seed = 12345
                else 
                    call system_clock(count=seed)
                end if 
                call random_seed(size=n)
                allocate(seed_array(n))
                seed_array = seed + 7*[(i-1,i=1,n)]
                call random_seed(put=seed_array)
                deallocate(seed_array)

                select case(trim(this%turb_intensity))
                case('light')
                    allocate(this%turb_hag(3))
                    allocate(this%turb_sig(3))
                    this%turb_hag = this%light_hag
                    this%turb_sig = this%light_sig 
                case('moderate')
                    allocate(this%turb_hag(3))
                    allocate(this%turb_sig(3))
                    this%turb_hag = this%moderate_hag
                    this%turb_sig = this%moderate_sig 
                case('severe')
                    allocate(this%turb_hag(4))
                    allocate(this%turb_sig(4))
                    this%turb_hag = this%severe_hag
                    this%turb_sig = this%severe_sig             
                end select 

                select case(trim(this%turb_model))
                    case('dryden_beal')
                        this%Lu = 1750.0
                        this%Lv = 875.0
                        this%Lw = 875.0
                        this%Lb = 4*this%wingspan/PI
                    case('dryden_8785')
                        this%Lu = 1750.0
                        this%Lv = 875.0
                        this%Lw = 875.0
                        this%Lb = 4*this%wingspan/PI
                end select 
                this%xff = 0.0
                this%prev_turb(:) = 0.0
                this%prev_xyz(:) = 0.0
                this%prev_f = 0.0
                this%prev_g = 0.0
                call json_get(j_turb, 'sample', j_sample, found)
                if(found) call turbulence_sample(this, j_sample)
            end if 
        end if 
    end subroutine atmosphere_init

    subroutine turbulence_sample(this, j_sample)
        implicit none 
        type(atmosphere_t) :: this
        type(json_value), pointer :: j_sample 
        character(len=:), allocatable :: fn 
        integer :: i, j, n, n_psd, iunit, psd_mean_unit 
        real, allocatable :: vals(:,:), psd_mean(:,:),psd_temp(:,:)
        real :: hag, sigma, dx, turb(6) ! hag means height above ground 
        real :: mean, stdev 
        logical :: found 
        ! Test random number generator function 
        ! call test_rand_normal()
        write(*,*) '    Sampling Atmospheric Turbulence...'
        call jsonx_get(j_sample, 'save_filename', fn)
        open(newunit=iunit, file=fn,status='REPLACE')
        write(*,*) '    saving sample to ', fn 
        call jsonx_get(j_sample, 'number_of_points', n) 
        call jsonx_get(j_sample, 'dx[ft]', dx)
        call jsonx_get(j_sample, 'height_above_ground[ft]', hag)
        allocate(vals(n,4))

        sigma = interpolate_1D(this%turb_hag, this%turb_sig, hag)
        write(*,*) '    Altitude = ', hag 
        write(*,*) '    Turbulence Standard Deviation, sigma = ', sigma
        write(iunit,*) 'distance[ft],uprime[ft/s],vprime[ft/s],wprime[ft/s],pprime[rad/s]' 
        do i = 1,n
            turb(:) = get_turbulence(this,dx,sigma,sigma,sigma)
            write(iunit,*) dx*real(i-1),',',turb(1),',',turb(2),',',turb(3),',',turb(4)
            vals(i,:) = turb(:)
        end do 
        close(iunit)

        call json_get(j_sample, 'psd_analyses', n_psd, found)
        if(found) then 
            call jsonx_get(j_sample, 'psd_analyses', n_psd)
            allocate(psd_mean(n/2+1,2))
            allocate(psd_temp(n/2+1,2))
            psd_mean = 0.0
            psd_temp = 0.0

            write(*,*) 'Turbulence Normalized Mean PSD Analysis'
            open(newunit=psd_mean_unit, file = 'PSD_Mean_Analysis.csv', status = 'REPLACE')
            write(*,*) '    - saving normalized PSD analysis to PSD_Mean_Analysis.csv'
            do j = 1, n_psd 
                write(*,*) 'PSD ',j,' of ',n_psd 
                do i =1,n
                    turb(:) = get_turbulence(this,dx,sigma,sigma,sigma)
                    vals(i,:) = turb(:) 
                end do 
                call psd(vals(:,4),dx,psd_norm=psd_temp) ! 1=u, 2=v, 3=w, 4=p
                psd_mean(:,2) = psd_mean(:,2) + psd_temp(:,2)/n_psd 
            end do 

            do i = 1,n/2+1
                write(psd_mean_unit,*) psd_temp(i,1)*2*PI,',',psd_mean(i,2)/2/PI ! x[rad/ft], y[ft]
            end do 
            close(psd_mean_unit)
        end if 
    end subroutine turbulence_sample

    function atmosphere_get_turbulence(this, states) result(ans)
        implicit none 
        type(atmosphere_t) :: this
        real :: states(21)
        real :: ans(4)
        real :: dx, sigma 
        dx = sqrt((states(7)-this%prev_xyz(1))**2 + (states(8)-this%prev_xyz(2))**2 + (states(9)-this%prev_xyz(3))**2)
        sigma = interpolate_1D(this%turb_hag, this%turb_sig, -states(9)) ! Assumes ground height is sea-level
        ans(:) = get_turbulence(this, dx, sigma, sigma, sigma)
        this%prev_xyz(:) = states(7:9)
    end function atmosphere_get_turbulence

    function get_turbulence(this,dx,su,sv,sw) result(ans)
        implicit none 
        type(atmosphere_t) :: this
        real :: dx, su, sv, sw 
        real :: ans(6)
        this%xff = this%xff + dx 
        select case(trim(this%turb_model))
            case('dryden_beal')
                ans(:) = dryden_beal(this,dx,su,sv,sw)
            case('dryden_8785')
                ans(:) = dryden_beal(this,dx,su,sv,sw) ! placeholder
        end select 
    end function get_turbulence

    function dryden_beal(this,dx,su,sv,sw) result(ans)
        implicit none 
        type(atmosphere_t) :: this
        real :: dx, su, sv, sw 
        real :: ans(6)
        real :: Au, Av, Aw, Ap 
        real :: etau, etav, etaw, etap 
        real :: f, g 

        Au = 0.5*dx/this%Lu 
        Av = 0.25*dx/this%Lv 
        Aw = 0.25*dx/this%Lw
        Ap = 0.5*dx/this%Lb

        etau = rand_normal()*su*sqrt(2.0*this%Lu/dx)
        etav = rand_normal()*sv*sqrt(2.0*this%Lv/dx)
        etaw = rand_normal()*sw*sqrt(2.0*this%Lw/dx)
        etap = rand_normal()*sw*sqrt(0.8*PI*(this%Lw/this%Lb)**(1.0/3.0)/this%Lw/dx)

        f = ((1.0-Av)*this%prev_f + 2.0*Av*etav)/(1.0+Av)
        g = ((1.0-Aw)*this%prev_g + 2.0*Aw*etaw)/(1.0+Aw)

        ans(1) = ((1.0-Au)*this%prev_turb(1)+2.0*Au*etau)/(1.0 + Au)
        ans(2) = ((1.0-Av)*this%prev_turb(2)+Av*(f+this%prev_f)+sqrt(3.)*(f-this%prev_f))/(1.0+Av)
        ans(3) = ((1.0-Aw)*this%prev_turb(3)+Aw*(g+this%prev_g)+sqrt(3.)*(g-this%prev_g))/(1.0+Aw)
        ans(4) = ((1.0-Ap)*this%prev_turb(4) + 2.0*Ap*etap)/(1 + Ap)
        ! now do disturbances in p,q,r here (ans(4), ans(5), ans(6))

        this%prev_f = f 
        this%prev_g = g 
        this%prev_turb(:) = ans(:)
    end function dryden_beal

    function dryden_8785(this,dx,su,sv,sw) result(ans) !!!! PLACEHOLDER FOR NOW !!!!
        implicit none 
        type(atmosphere_t) :: this
        real :: dx, su, sv, sw 
        real :: ans(4)
        ! AS PLACEHOLDER, CALL BEAL VERSION UNTIL THIS IS IMPLEMENTED
        ans = dryden_beal(this,dx,su,sv,sw)
    end function dryden_8785

end module atmosphere_m