module control
    implicit none
    character   slash

    !da_get_unit
    integer,parameter:: unit_start=1000
    integer,parameter:: unit_end=3000
    logical:: unit_used(unit_start:unit_end)=.false.

    !namelist
    !&system_setting
    character*7 system_type !Choose from 'Windows' or 'Linux'
    integer max_openmp_threads
    namelist /system_setting/ system_type,max_openmp_threads

    real run_hours,history_interval
    namelist /time_control/ run_hours,history_interval

    !&files
    character*500 input_data_path     ,&
                  output_data_path

    namelist /files/ input_data_path    ,&
                     output_data_path

    !&domains
    integer::   s_we=1,&
                s_sn=1,&
                e_we  ,&
                e_sn

    !&dynamic
    integer ::  integration_option = 1
    real    ::  z0                 = 2500.0            ,&
                smooth_coefficient = 0.5               ,&
                spec_exp           = 0.3333333333333333
    namelist /dynamics/ integration_option,z0,smooth_coefficient,spec_exp

    real time_step

    namelist /domains/ time_step

    !output
    integer output_time_num
    real history_interval_seconds
    character*33 XTIME_units

    !integration
    type domain
        integer         x_grd_num,y_grd_num     ,&
                        x_grd_num_u,y_grd_num_u ,&
                        x_grd_num_v,y_grd_num_v

        integer         map_proj,bdy_width,bdy_time_num,bdy_tend_time_num

        real            integerated_step_num,interval_seconds,dx,cen_lat,cen_lon,truelat1,truelat2,stand_lon

        character,allocatable:: times(:,:),bdy_times(:,:)
        integer start_year      ,&
                start_month     ,&
                start_day       ,&
                start_hour      ,&
                start_minute    ,&
                start_second    ,&
                end_year        ,&
                end_month       ,&
                end_day         ,&
                end_hour        ,&
                end_minute      ,&
                end_second      ,&
                spec_zone       ,&
                relax_zone

        real,allocatable::      z(:,:,:)        ,&
                                u(:,:,:)        ,&
                                v(:,:,:)        ,&
                                ua(:,:)       ,&
                                ub(:,:)       ,&
                                uc(:,:)       ,&
                                va(:,:)       ,&
                                vb(:,:)       ,&
                                vc(:,:)       ,&
                                za(:,:)       ,&
                                zb(:,:)       ,&
                                zc(:,:)       ,&
                                mapfac_m(:,:)   ,&
                                mapfac_c(:,:)   ,&
                                mapfac_u(:,:)   ,&
                                mapfac_v(:,:)   ,&
                                f_m(:,:)        ,&
                                f_c(:,:)        ,&
                                f_u(:,:)        ,&
                                f_v(:,:)        ,&
                                XTIME(:)        ,&
                                XLAT_M(:,:)     ,&
                                XLONG_M(:,:)    ,&
                                XLAT_C(:,:)     ,&
                                XLONG_C(:,:)    ,&
                                XLAT_U(:,:)     ,&
                                XLONG_U(:,:)    ,&
                                XLAT_V(:,:)     ,&
                                XLONG_V(:,:)    ,&
                                Z_BYS(:,:,:)    ,&
                                Z_BYE(:,:,:)    ,&
                                Z_BXS(:,:,:)    ,&
                                Z_BXE(:,:,:)    ,&
                                U_BYS(:,:,:)    ,&
                                U_BYE(:,:,:)    ,&
                                U_BXS(:,:,:)    ,&
                                U_BXE(:,:,:)    ,&
                                V_BYS(:,:,:)    ,&
                                V_BYE(:,:,:)    ,&
                                V_BXS(:,:,:)    ,&
                                V_BXE(:,:,:)    ,&
                                Z_BTYS(:,:,:)   ,&
                                Z_BTYE(:,:,:)   ,&
                                Z_BTXS(:,:,:)   ,&
                                Z_BTXE(:,:,:)   ,&
                                U_BTYS(:,:,:)   ,&
                                U_BTYE(:,:,:)   ,&
                                U_BTXS(:,:,:)   ,&
                                U_BTXE(:,:,:)   ,&
                                V_BTYS(:,:,:)   ,&
                                V_BTYE(:,:,:)   ,&
                                V_BTXS(:,:,:)   ,&
                                V_BTXE(:,:,:)   ,&
                                F1(:)           ,&
                                F2(:)           ,&
                                K(:)
    endtype

    type(domain)::  grid
    
    !Boundary
    real::  filling_value=-999999999.99d0
    
    !Timer
    real :: time_divide_value

    contains
    subroutine da_get_unit(unit)
        implicit none

        integer, intent(out) :: unit

        integer :: i

        unit = -1

        do i = unit_start, unit_end
           if (.NOT. unit_used(i)) then
              unit=i
              unit_used(i) = .true.
              exit
           end if
        end do

        if (unit == -1) then
          write(*,*)"No free units"
        end if
    end subroutine da_get_unit

    subroutine read_namelist
        implicit none
        integer unit

        call da_get_unit(unit)
        open(unit,file='namelist.input',status='old')
        read(unit,nml=system_setting)
        read(unit,nml=time_control)
        read(unit,nml=files)
        read(unit,nml=domains)
        read(unit,nml=dynamics)
        close(unit)

    end subroutine read_namelist

    subroutine choose_slash
    !This subroutine must be used after using read_namelist
        implicit none
        if(trim(adjustl(system_type))=='Windows'.or.trim(adjustl(system_type))=='windows')then
            slash='\'
            time_divide_value=10000.
        elseif(trim(adjustl(system_type))=='Linux'.or.trim(adjustl(system_type))=='linux')then
            slash='/'
            time_divide_value=1000.
        endif
    end subroutine choose_slash

end module control