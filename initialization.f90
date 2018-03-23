module initialization
use control
use netcdf
use tools
implicit none
include 'netcdf.inc'

contains
    subroutine allocate_arrays
        allocate(   grid%z          (grid%x_grd_num     ,grid%y_grd_num     ,grid%BDY_TIME_NUM) ,&
                    grid%u          (grid%x_grd_num_u   ,grid%y_grd_num_u   ,grid%BDY_TIME_NUM) ,&
                    grid%v          (grid%x_grd_num_v   ,grid%y_grd_num_v   ,grid%BDY_TIME_NUM) ,&
                    grid%ua         (grid%x_grd_num_u   ,grid%y_grd_num_u)                      ,&
                    grid%ub         (grid%x_grd_num_u   ,grid%y_grd_num_u)                      ,&
                    grid%uc         (grid%x_grd_num_u   ,grid%y_grd_num_u)                      ,&
                    grid%va         (grid%x_grd_num_v   ,grid%y_grd_num_v)                      ,&
                    grid%vb         (grid%x_grd_num_v   ,grid%y_grd_num_v)                      ,&
                    grid%vc         (grid%x_grd_num_v   ,grid%y_grd_num_v)                      ,&
                    grid%za         (grid%x_grd_num     ,grid%y_grd_num)                        ,&
                    grid%zb         (grid%x_grd_num     ,grid%y_grd_num)                        ,&
                    grid%zc         (grid%x_grd_num     ,grid%y_grd_num)                        ,&
                    grid%mapfac_m   (grid%x_grd_num     ,grid%y_grd_num)                        ,&
                    grid%mapfac_c   (grid%x_grd_num_u   ,grid%y_grd_num_v)                      ,&
                    grid%mapfac_u   (grid%x_grd_num_u   ,grid%y_grd_num_u)                      ,&
                    grid%mapfac_v   (grid%x_grd_num_v   ,grid%y_grd_num_v)                      ,&
                    grid%f_m        (grid%x_grd_num     ,grid%y_grd_num)                        ,&
                    grid%f_c        (grid%x_grd_num_u   ,grid%y_grd_num_v)                      ,&
                    grid%f_u        (grid%x_grd_num_u   ,grid%y_grd_num_u)                      ,&
                    grid%f_v        (grid%x_grd_num_v   ,grid%y_grd_num_v)                      ,&
                    grid%XLAT_M     (grid%x_grd_num     ,grid%y_grd_num)                        ,&
                    grid%XLONG_M    (grid%x_grd_num     ,grid%y_grd_num)                        ,&
                    grid%XLAT_C     (grid%x_grd_num_u   ,grid%y_grd_num_v)                      ,&
                    grid%XLONG_C    (grid%x_grd_num_u   ,grid%y_grd_num_v)                      ,&
                    grid%XLAT_U     (grid%x_grd_num_u   ,grid%y_grd_num_u)                      ,&
                    grid%XLONG_U    (grid%x_grd_num_u   ,grid%y_grd_num_u)                      ,&
                    grid%XLAT_V     (grid%x_grd_num_v   ,grid%y_grd_num_v)                      ,&
                    grid%XLONG_V    (grid%x_grd_num_v   ,grid%y_grd_num_v)                      ,&

                    grid%Z_BYS   (grid%x_grd_num    ,grid%bdy_width     ,grid%BDY_TIME_NUM)     ,&
                    grid%Z_BYE   (grid%x_grd_num    ,grid%bdy_width     ,grid%BDY_TIME_NUM)     ,&
                    grid%Z_BXS   (grid%y_grd_num    ,grid%bdy_width     ,grid%BDY_TIME_NUM)     ,&
                    grid%Z_BXE   (grid%y_grd_num    ,grid%bdy_width     ,grid%BDY_TIME_NUM)     ,&

                    grid%U_BYS   (grid%x_grd_num_u  ,grid%bdy_width     ,grid%BDY_TIME_NUM)     ,&
                    grid%U_BYE   (grid%x_grd_num_u  ,grid%bdy_width     ,grid%BDY_TIME_NUM)     ,&
                    grid%U_BXS   (grid%y_grd_num_u  ,grid%bdy_width     ,grid%BDY_TIME_NUM)     ,&
                    grid%U_BXE   (grid%y_grd_num_u  ,grid%bdy_width     ,grid%BDY_TIME_NUM)     ,&

                    grid%V_BYS   (grid%x_grd_num_v  ,grid%bdy_width     ,grid%BDY_TIME_NUM)     ,&
                    grid%V_BYE   (grid%x_grd_num_v  ,grid%bdy_width     ,grid%BDY_TIME_NUM)     ,&
                    grid%V_BXS   (grid%y_grd_num_v  ,grid%bdy_width     ,grid%BDY_TIME_NUM)     ,&
                    grid%V_BXE   (grid%y_grd_num_v  ,grid%bdy_width     ,grid%BDY_TIME_NUM)     ,&

                    grid%Z_BTYS  (grid%x_grd_num    ,grid%bdy_width     ,grid%BDY_TEND_TIME_NUM),&
                    grid%Z_BTYE  (grid%x_grd_num    ,grid%bdy_width     ,grid%BDY_TEND_TIME_NUM),&
                    grid%Z_BTXS  (grid%y_grd_num    ,grid%bdy_width     ,grid%BDY_TEND_TIME_NUM),&
                    grid%Z_BTXE  (grid%y_grd_num    ,grid%bdy_width     ,grid%BDY_TEND_TIME_NUM),&

                    grid%U_BTYS  (grid%x_grd_num_u  ,grid%bdy_width     ,grid%BDY_TEND_TIME_NUM),&
                    grid%U_BTYE  (grid%x_grd_num_u  ,grid%bdy_width     ,grid%BDY_TEND_TIME_NUM),&
                    grid%U_BTXS  (grid%y_grd_num_u  ,grid%bdy_width     ,grid%BDY_TEND_TIME_NUM),&
                    grid%U_BTXE  (grid%y_grd_num_u  ,grid%bdy_width     ,grid%BDY_TEND_TIME_NUM),&

                    grid%V_BTYS  (grid%x_grd_num_v  ,grid%bdy_width     ,grid%BDY_TEND_TIME_NUM),&
                    grid%V_BTYE  (grid%x_grd_num_v  ,grid%bdy_width     ,grid%BDY_TEND_TIME_NUM),&
                    grid%V_BTXS  (grid%y_grd_num_v  ,grid%bdy_width     ,grid%BDY_TEND_TIME_NUM),&
                    grid%V_BTXE  (grid%y_grd_num_v  ,grid%bdy_width     ,grid%BDY_TEND_TIME_NUM),&

                    grid%XTIME    (output_time_num      )                                       ,&
                    grid%times    (19,output_time_num   )                                       ,&
                    grid%bdy_times(19,grid%BDY_TIME_NUM )                                       ,&
                    grid%F1       (grid%Spec_Zone+1:grid%bdy_width-1)                           ,&
                    grid%F2       (grid%Spec_Zone+1:grid%bdy_width-1)                           ,&
                    grid%K        (1:grid%bdy_width)                                            )
                    
    end subroutine allocate_arrays

    subroutine init_data
        implicit none
        !netCDF
        integer ncid,status

        !Variables ID
        integer bdy_times_id,&
                XLAT_M_id   ,&
                XLONG_M_id  ,&
                XLAT_C_id   ,&
                XLONG_C_id  ,&
                XLAT_U_id   ,&
                XLONG_U_id  ,&
                XLAT_V_id   ,&
                XLONG_V_id  ,&
                MAPFAC_M_id ,&
                MAPFAC_C_id ,&
                MAPFAC_U_id ,&
                MAPFAC_V_id ,&
                F_M_id      ,&
                F_C_id      ,&
                F_U_id      ,&
                F_V_id      ,&
                Z_id        ,&
                U_id        ,&
                V_id        ,&
                Z_BYS_id    ,&
                Z_BYE_id    ,&
                Z_BXS_id    ,&
                Z_BXE_id    ,&
                U_BYS_id    ,&
                U_BYE_id    ,&
                U_BXS_id    ,&
                U_BXE_id    ,&
                V_BYS_id    ,&
                V_BYE_id    ,&
                V_BXS_id    ,&
                V_BXE_id    ,&
                Z_BTYS_id   ,&
                Z_BTYE_id   ,&
                Z_BTXS_id   ,&
                Z_BTXE_id   ,&
                U_BTYS_id   ,&
                U_BTYE_id   ,&
                U_BTXS_id   ,&
                U_BTXE_id   ,&
                V_BTYS_id   ,&
                V_BTYE_id   ,&
                V_BTXS_id   ,&
                V_BTXE_id

        status=nf90_open(trim(adjustl(input_data_path))//slash//'BPEM_input.nc',NF_NOWRITE,ncid)
        status=nf90_get_att(ncid,NF90_GLOBAL,'MAP_PROJ'                  ,grid%map_proj         )
        status=nf90_get_att(ncid,NF90_GLOBAL,'BDY_WIDTH'                 ,grid%bdy_width        )
        status=nf90_get_att(ncid,NF90_GLOBAL,'BDY_TIME_NUM'              ,grid%bdy_time_num     )
        status=nf90_get_att(ncid,NF90_GLOBAL,'BDY_TEND_TIME_NUM'         ,grid%bdy_tend_time_num)
        status=nf90_get_att(ncid,NF90_GLOBAL,'WEST-EAST_GRID_DIMENSION'  ,grid%x_grd_num_u      )
        status=nf90_get_att(ncid,NF90_GLOBAL,'SOUTH-NORTH_GRID_DIMENSION',grid%y_grd_num_v      )
        status=nf90_get_att(ncid,NF90_GLOBAL,'DX'                        ,grid%dx               )
        status=nf90_get_att(ncid,NF90_GLOBAL,'CEN_LAT'                   ,grid%CEN_LAT          )
        status=nf90_get_att(ncid,NF90_GLOBAL,'CEN_LON'                   ,grid%CEN_LON          )
        status=nf90_get_att(ncid,NF90_GLOBAL,'TRUELAT1'                  ,grid%TRUELAT1         )
        status=nf90_get_att(ncid,NF90_GLOBAL,'TRUELAT2'                  ,grid%TRUELAT2         )
        status=nf90_get_att(ncid,NF90_GLOBAL,'STAND_LON'                 ,grid%STAND_LON        )
        status=nf90_get_att(ncid,NF90_GLOBAL,'start_year'                ,grid%start_year       )
        status=nf90_get_att(ncid,NF90_GLOBAL,'start_month'               ,grid%start_month      )
        status=nf90_get_att(ncid,NF90_GLOBAL,'start_day'                 ,grid%start_day        )
        status=nf90_get_att(ncid,NF90_GLOBAL,'start_hour'                ,grid%start_hour       )
        status=nf90_get_att(ncid,NF90_GLOBAL,'start_minute'              ,grid%start_minute     )
        status=nf90_get_att(ncid,NF90_GLOBAL,'start_second'              ,grid%start_second     )
        status=nf90_get_att(ncid,NF90_GLOBAL,'end_year'                  ,grid%end_year         )
        status=nf90_get_att(ncid,NF90_GLOBAL,'end_month'                 ,grid%end_month        )
        status=nf90_get_att(ncid,NF90_GLOBAL,'end_day'                   ,grid%end_day          )
        status=nf90_get_att(ncid,NF90_GLOBAL,'end_hour'                  ,grid%end_hour         )
        status=nf90_get_att(ncid,NF90_GLOBAL,'end_minute'                ,grid%end_minute       )
        status=nf90_get_att(ncid,NF90_GLOBAL,'end_second'                ,grid%end_second       )
        status=nf90_get_att(ncid,NF90_GLOBAL,'interval_seconds'          ,grid%interval_seconds )
        status=nf90_get_att(ncid,NF90_GLOBAL,'spec_zone'                 ,grid%spec_zone        )
        status=nf90_get_att(ncid,NF90_GLOBAL,'relax_zone'                ,grid%relax_zone       )

        grid%x_grd_num  =grid%x_grd_num_u-1
        grid%y_grd_num  =grid%y_grd_num_v-1
        grid%y_grd_num_u=grid%y_grd_num
        grid%x_grd_num_v=grid%x_grd_num

        output_time_num=ceiling(run_hours*60/history_interval)+1
        history_interval_seconds=history_interval*60
        call allocate_arrays

        status=nf90_inq_varid(ncid,'bdy_times',bdy_times_id)
        status=nf90_inq_varid(ncid,'XLAT_M'   ,XLAT_M_id   )
        status=nf90_inq_varid(ncid,'XLONG_M'  ,XLONG_M_id  )
        status=nf90_inq_varid(ncid,'XLAT_C'   ,XLAT_C_id   )
        status=nf90_inq_varid(ncid,'XLONG_C'  ,XLONG_C_id  )
        status=nf90_inq_varid(ncid,'XLAT_U'   ,XLAT_U_id   )
        status=nf90_inq_varid(ncid,'XLONG_U'  ,XLONG_U_id  )
        status=nf90_inq_varid(ncid,'XLAT_V'   ,XLAT_V_id   )
        status=nf90_inq_varid(ncid,'XLONG_V'  ,XLONG_V_id  )
        status=nf90_inq_varid(ncid,'MAPFAC_M' ,MAPFAC_M_id )
        status=nf90_inq_varid(ncid,'MAPFAC_C' ,MAPFAC_C_id )
        status=nf90_inq_varid(ncid,'MAPFAC_U' ,MAPFAC_U_id )
        status=nf90_inq_varid(ncid,'MAPFAC_V' ,MAPFAC_V_id )
        status=nf90_inq_varid(ncid,'F_M'      ,F_M_id      )
        status=nf90_inq_varid(ncid,'F_C'      ,F_C_id      )
        status=nf90_inq_varid(ncid,'F_U'      ,F_U_id      )
        status=nf90_inq_varid(ncid,'F_V'      ,F_V_id      )
        status=nf90_inq_varid(ncid,'Z'        ,Z_id        )
        status=nf90_inq_varid(ncid,'U'        ,U_id        )
        status=nf90_inq_varid(ncid,'V'        ,V_id        )
        status=nf90_inq_varid(ncid,'Z_BYS'    ,Z_BYS_id    )
        status=nf90_inq_varid(ncid,'Z_BYE'    ,Z_BYE_id    )
        status=nf90_inq_varid(ncid,'Z_BXS'    ,Z_BXS_id    )
        status=nf90_inq_varid(ncid,'Z_BXE'    ,Z_BXE_id    )
        status=nf90_inq_varid(ncid,'U_BYS'    ,U_BYS_id    )
        status=nf90_inq_varid(ncid,'U_BYE'    ,U_BYE_id    )
        status=nf90_inq_varid(ncid,'U_BXS'    ,U_BXS_id    )
        status=nf90_inq_varid(ncid,'U_BXE'    ,U_BXE_id    )
        status=nf90_inq_varid(ncid,'V_BYS'    ,V_BYS_id    )
        status=nf90_inq_varid(ncid,'V_BYE'    ,V_BYE_id    )
        status=nf90_inq_varid(ncid,'V_BXS'    ,V_BXS_id    )
        status=nf90_inq_varid(ncid,'V_BXE'    ,V_BXE_id    )
        status=nf90_inq_varid(ncid,'Z_BTYS'   ,Z_BTYS_id   )
        status=nf90_inq_varid(ncid,'Z_BTYE'   ,Z_BTYE_id   )
        status=nf90_inq_varid(ncid,'Z_BTXS'   ,Z_BTXS_id   )
        status=nf90_inq_varid(ncid,'Z_BTXE'   ,Z_BTXE_id   )
        status=nf90_inq_varid(ncid,'U_BTYS'   ,U_BTYS_id   )
        status=nf90_inq_varid(ncid,'U_BTYE'   ,U_BTYE_id   )
        status=nf90_inq_varid(ncid,'U_BTXS'   ,U_BTXS_id   )
        status=nf90_inq_varid(ncid,'U_BTXE'   ,U_BTXE_id   )
        status=nf90_inq_varid(ncid,'V_BTYS'   ,V_BTYS_id   )
        status=nf90_inq_varid(ncid,'V_BTYE'   ,V_BTYE_id   )
        status=nf90_inq_varid(ncid,'V_BTXS'   ,V_BTXS_id   )
        status=nf90_inq_varid(ncid,'V_BTXE'   ,V_BTXE_id   )

        status=nf90_get_var(ncid,bdy_times_id,grid%bdy_times,start=(/1,1/),count=(/19,grid%bdy_time_num/))
        status=nf90_get_var(ncid,XLAT_M_id   ,grid%XLAT_M    )
        status=nf90_get_var(ncid,XLONG_M_id  ,grid%XLONG_M   )
        status=nf90_get_var(ncid,XLAT_C_id   ,grid%XLAT_C    )
        status=nf90_get_var(ncid,XLONG_C_id  ,grid%XLONG_C   )
        status=nf90_get_var(ncid,XLAT_U_id   ,grid%XLAT_U    )
        status=nf90_get_var(ncid,XLONG_U_id  ,grid%XLONG_U   )
        status=nf90_get_var(ncid,XLAT_V_id   ,grid%XLAT_V    )
        status=nf90_get_var(ncid,XLONG_V_id  ,grid%XLONG_V   )
        status=nf90_get_var(ncid,MAPFAC_M_id ,grid%MAPFAC_M  )
        status=nf90_get_var(ncid,MAPFAC_C_id ,grid%MAPFAC_C  )
        status=nf90_get_var(ncid,MAPFAC_U_id ,grid%MAPFAC_U  )
        status=nf90_get_var(ncid,MAPFAC_V_id ,grid%MAPFAC_V  )
        status=nf90_get_var(ncid,F_M_id      ,grid%F_M       )
        status=nf90_get_var(ncid,F_C_id      ,grid%F_C       )
        status=nf90_get_var(ncid,F_U_id      ,grid%F_U       )
        status=nf90_get_var(ncid,F_V_id      ,grid%F_V       )
        status=nf90_get_var(ncid,Z_id        ,grid%Z         )
        status=nf90_get_var(ncid,U_id        ,grid%U         )
        status=nf90_get_var(ncid,V_id        ,grid%V         )
        status=nf90_get_var(ncid,Z_BYS_id    ,grid%Z_BYS     )
        status=nf90_get_var(ncid,Z_BYE_id    ,grid%Z_BYE     )
        status=nf90_get_var(ncid,Z_BXS_id    ,grid%Z_BXS     )
        status=nf90_get_var(ncid,Z_BXE_id    ,grid%Z_BXE     )
        status=nf90_get_var(ncid,U_BYS_id    ,grid%U_BYS     )
        status=nf90_get_var(ncid,U_BYE_id    ,grid%U_BYE     )
        status=nf90_get_var(ncid,U_BXS_id    ,grid%U_BXS     )
        status=nf90_get_var(ncid,U_BXE_id    ,grid%U_BXE     )
        status=nf90_get_var(ncid,V_BYS_id    ,grid%V_BYS     )
        status=nf90_get_var(ncid,V_BYE_id    ,grid%V_BYE     )
        status=nf90_get_var(ncid,V_BXS_id    ,grid%V_BXS     )
        status=nf90_get_var(ncid,V_BXE_id    ,grid%V_BXE     )
        status=nf90_get_var(ncid,Z_BTYS_id   ,grid%Z_BTYS    )
        status=nf90_get_var(ncid,Z_BTYE_id   ,grid%Z_BTYE    )
        status=nf90_get_var(ncid,Z_BTXS_id   ,grid%Z_BTXS    )
        status=nf90_get_var(ncid,Z_BTXE_id   ,grid%Z_BTXE    )
        status=nf90_get_var(ncid,U_BTYS_id   ,grid%U_BTYS    )
        status=nf90_get_var(ncid,U_BTYE_id   ,grid%U_BTYE    )
        status=nf90_get_var(ncid,U_BTXS_id   ,grid%U_BTXS    )
        status=nf90_get_var(ncid,U_BTXE_id   ,grid%U_BTXE    )
        status=nf90_get_var(ncid,V_BTYS_id   ,grid%V_BTYS    )
        status=nf90_get_var(ncid,V_BTYE_id   ,grid%V_BTYE    )
        status=nf90_get_var(ncid,V_BTXS_id   ,grid%V_BTXS    )
        status=nf90_get_var(ncid,V_BTXE_id   ,grid%V_BTXE    )

        status=nf_close(ncid)

        grid%za=grid%z(:,:,1)
        grid%ua=grid%u(:,:,1)
        grid%va=grid%v(:,:,1)
        grid%zb=0
        grid%ub=0
        grid%vb=0
        grid%zc=0
        grid%uc=0
        grid%vc=0

        call calculate_output_times

    end subroutine init_data

    !Calculate output times for grid%times, this program has to be run after knowing the output_time_num
    subroutine calculate_output_times
    implicit none
    integer year,month,day,hour,minute,second
    character*4 year_char
    character*2 month_char,day_char,hour_char,minute_char,second_char
    character*16 time_increase_option
    integer time_increase_value
    character*19 times_temp
    integer i,j

    year    =   grid%start_year
    month   =   grid%start_month
    day     =   grid%start_day
    hour    =   grid%start_hour
    minute  =   grid%start_minute
    second  =   grid%start_second

    time_increase_option='calculate_minute'
    time_increase_value=int(history_interval)
    do i=1,output_time_num
        write(year_char     ,'(i4.4)')year
        write(month_char    ,'(i2.2)')month
        write(day_char      ,'(i2.2)')day
        write(hour_char     ,'(i2.2)')hour
        write(minute_char   ,'(i2.2)')minute
        write(second_char   ,'(i2.2)')second

        if(i==1)XTIME_units='minutes since '//year_char//'-'//month_char//'-'//day_char//' '//hour_char//':'//minute_char//':'//second_char

        times_temp=year_char//'-'//month_char//'-'//day_char//'_'//hour_char//':'//minute_char//':'//second_char

        do j=1,19
            grid%times(j,i)=times_temp(j:j)
        enddo

        call time_calculator(year,month,day,hour,minute,second,&
                             year,month,day,hour,minute,second,&
                             time_increase_option,time_increase_value)

        grid%XTIME(i)=(i-1)*history_interval
    enddo


    end subroutine calculate_output_times

end module initialization