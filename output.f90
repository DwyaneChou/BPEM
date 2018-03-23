module output
    use control
    use netcdf
    use tools
    implicit none
    include 'netcdf.inc'
    contains

    !Create netCDF output file
    subroutine create_output
    implicit none
    integer ncid,status
    !Dimensions ID
    integer west_east_dimid         ,&
            west_east_stag_dimid    ,&
            south_north_dimid       ,&
            south_north_stag_dimid  ,&
            Time_dimid              ,&
            DateStrLen_dimid

    !Variables ID
    integer Time_id    ,&
            Times_id   ,&
            XTIME_id   ,&
            XLAT_M_id  ,&
            XLONG_M_id ,&
            XLAT_U_id  ,&
            XLONG_U_id ,&
            XLAT_V_id  ,&
            XLONG_V_id ,&
            U_id       ,&
            V_id       ,&
            Z_id

    status=NF90_CREATE(trim(adjustl(output_data_path))//slash//'BPEM_output.nc',OR(NF_SHARE,NF_64BIT_OFFSET),ncid)

    !Get dimension IDs
    status=NF90_DEF_DIM(ncid,'west_east'        ,grid%x_grd_num     ,west_east_dimid        )
    status=NF90_DEF_DIM(ncid,'west_east_stag'   ,grid%x_grd_num_u   ,west_east_stag_dimid   )
    status=NF90_DEF_DIM(ncid,'south_north'      ,grid%y_grd_num     ,south_north_dimid      )
    status=NF90_DEF_DIM(ncid,'south_north_stag' ,grid%y_grd_num_v   ,south_north_stag_dimid )
    status=NF90_DEF_DIM(ncid,'Time'             ,output_time_num    ,Time_dimid             )
    status=NF90_DEF_DIM(ncid,'DateStrLen'       ,19                 ,DateStrLen_dimid       )

    !Define Variables
    status=NF90_DEF_VAR(ncid,'Time'     ,NF90_DOUBLE,(/Time_dimid/)                                         ,Time_id   )
    status=NF90_PUT_ATT(ncid,Time_id    ,'units'       ,XTIME_units)
    status=NF90_PUT_ATT(ncid,Time_id    ,'description' ,XTIME_units)

    status=NF90_DEF_VAR(ncid,'Times'    ,NF90_CHAR  ,(/DateStrLen_dimid,Time_dimid/)                        ,Times_id   )

    status=NF90_DEF_VAR(ncid,'XTIME'    ,NF90_DOUBLE,(/Time_dimid/)                                         ,XTIME_id   )
    status=NF90_PUT_ATT(ncid,XTIME_id,'units'       ,XTIME_units)
    status=NF90_PUT_ATT(ncid,XTIME_id,'description' ,XTIME_units)

    status=NF90_DEF_VAR(ncid,'XLAT_M'   ,NF90_DOUBLE,(/west_east_dimid,south_north_dimid,Time_dimid/)       ,XLAT_M_id  )
    status=NF90_PUT_ATT(ncid,XLAT_M_id,'units','degree_north')
    status=NF90_PUT_ATT(ncid,XLAT_M_id,'description','Latitude on mass grid')
    status=NF90_PUT_ATT(ncid,XLAT_M_id,'coordinates','XLONG_M XLAT_M')

    status=NF90_DEF_VAR(ncid,'XLONG_M'  ,NF90_DOUBLE,(/west_east_dimid,south_north_dimid,Time_dimid/)       ,XLONG_M_id )
    status=NF90_PUT_ATT(ncid,XLONG_M_id,'units','degree_east')
    status=NF90_PUT_ATT(ncid,XLONG_M_id,'description','Longitude on mass grid')
    status=NF90_PUT_ATT(ncid,XLONG_M_id,'coordinates','XLONG_M XLAT_M')

    status=NF90_DEF_VAR(ncid,'XLAT_U'   ,NF90_DOUBLE,(/west_east_stag_dimid,south_north_dimid,Time_dimid/)  ,XLAT_U_id  )
    status=NF90_PUT_ATT(ncid,XLAT_U_id,'units','degree_north')
    status=NF90_PUT_ATT(ncid,XLAT_U_id,'description','Latitude on U grid')
    status=NF90_PUT_ATT(ncid,XLAT_U_id,'coordinates','XLONG_U XLAT_U')

    status=NF90_DEF_VAR(ncid,'XLONG_U'  ,NF90_DOUBLE,(/west_east_stag_dimid,south_north_dimid,Time_dimid/)  ,XLONG_U_id )
    status=NF90_PUT_ATT(ncid,XLONG_U_id,'units','degree_east')
    status=NF90_PUT_ATT(ncid,XLONG_U_id,'description','Longitude on U grid')
    status=NF90_PUT_ATT(ncid,XLONG_U_id,'coordinates','XLONG_U XLAT_U')

    status=NF90_DEF_VAR(ncid,'XLAT_V'   ,NF90_DOUBLE,(/west_east_dimid,south_north_stag_dimid,Time_dimid/)  ,XLAT_V_id  )
    status=NF90_PUT_ATT(ncid,XLAT_V_id,'units','degree_north')
    status=NF90_PUT_ATT(ncid,XLAT_V_id,'description','Latitude on V grid')
    status=NF90_PUT_ATT(ncid,XLAT_V_id,'coordinates','XLONG_V XLAT_V')

    status=NF90_DEF_VAR(ncid,'XLONG_V'  ,NF90_DOUBLE,(/west_east_dimid,south_north_stag_dimid,Time_dimid/)  ,XLONG_V_id )
    status=NF90_PUT_ATT(ncid,XLONG_V_id,'units','degree_east')
    status=NF90_PUT_ATT(ncid,XLONG_V_id,'description','Longitude on V grid')
    status=NF90_PUT_ATT(ncid,XLONG_V_id,'coordinates','XLONG_V XLAT_V')

    status=NF90_DEF_VAR(ncid,'U'        ,NF90_DOUBLE,(/west_east_stag_dimid,south_north_dimid,Time_dimid/)  ,U_id       )
    status=NF90_PUT_ATT(ncid,U_id,'units','m/s')
    status=NF90_PUT_ATT(ncid,U_id,'description','U-Component Wind')
    status=NF90_PUT_ATT(ncid,U_id,'coordinates','XLONG_U XLAT_U XTIME')

    status=NF90_DEF_VAR(ncid,'V'        ,NF90_DOUBLE,(/west_east_dimid,south_north_stag_dimid,Time_dimid/)  ,V_id       )
    status=NF90_PUT_ATT(ncid,V_id,'units','m/s')
    status=NF90_PUT_ATT(ncid,V_id,'description','V-Component Wind')
    status=NF90_PUT_ATT(ncid,V_id,'coordinates','XLONG_V XLAT_V XTIME')

    status=NF90_DEF_VAR(ncid,'Z'        ,NF90_DOUBLE,(/west_east_dimid,south_north_dimid,Time_dimid/)       ,Z_id       )
    status=NF90_PUT_ATT(ncid,Z_id,'units','dagpm')
    status=NF90_PUT_ATT(ncid,Z_id,'description','Geopotential Height')
    status=NF90_PUT_ATT(ncid,Z_id,'coordinates','XLONG_M XLAT_M XTIME')

    !Put Global Attributes
    status=nf90_put_att(ncid,NF90_GLOBAL,'MAP_PROJ'                  ,grid%map_proj    )
    status=nf90_put_att(ncid,NF90_GLOBAL,'WEST-EAST_GRID_DIMENSION'  ,grid%x_grd_num_u )
    status=nf90_put_att(ncid,NF90_GLOBAL,'SOUTH-NORTH_GRID_DIMENSION',grid%y_grd_num_v )
    status=nf90_put_att(ncid,NF90_GLOBAL,'DX'                        ,grid%dx          )
    status=nf90_put_att(ncid,NF90_GLOBAL,'DY'                        ,grid%dx          )
    status=nf90_put_att(ncid,NF90_GLOBAL,'CEN_LAT'                   ,grid%CEN_LAT     )
    status=nf90_put_att(ncid,NF90_GLOBAL,'CEN_LON'                   ,grid%CEN_LON     )
    status=nf90_put_att(ncid,NF90_GLOBAL,'TRUELAT1'                  ,grid%TRUELAT1    )
    status=nf90_put_att(ncid,NF90_GLOBAL,'TRUELAT2'                  ,grid%TRUELAT2    )
    status=nf90_put_att(ncid,NF90_GLOBAL,'STAND_LON'                 ,grid%STAND_LON   )
    status=nf90_put_att(ncid,NF90_GLOBAL,'start_year'                ,grid%start_year  )
    status=nf90_put_att(ncid,NF90_GLOBAL,'start_month'               ,grid%start_month )
    status=nf90_put_att(ncid,NF90_GLOBAL,'start_day'                 ,grid%start_day   )
    status=nf90_put_att(ncid,NF90_GLOBAL,'start_hour'                ,grid%start_hour  )
    status=nf90_put_att(ncid,NF90_GLOBAL,'start_minute'              ,grid%start_minute)
    status=nf90_put_att(ncid,NF90_GLOBAL,'start_second'              ,grid%start_second)
    status=nf90_put_att(ncid,NF90_GLOBAL,'end_year'                  ,grid%end_year    )
    status=nf90_put_att(ncid,NF90_GLOBAL,'end_month'                 ,grid%end_month   )
    status=nf90_put_att(ncid,NF90_GLOBAL,'end_day'                   ,grid%end_day     )
    status=nf90_put_att(ncid,NF90_GLOBAL,'end_hour'                  ,grid%end_hour    )
    status=nf90_put_att(ncid,NF90_GLOBAL,'end_minute'                ,grid%end_minute  )
    status=nf90_put_att(ncid,NF90_GLOBAL,'end_second'                ,grid%end_second  )

    status=NF90_ENDDEF(ncid)

    !Put variables' value into the file
    print*,'It is time to write the output at ',grid%times(:,1)
    status=NF90_PUT_VAR(ncid,Times_id   ,grid%times,start=(/1,1/),count=(/19,output_time_num/))
    status=NF90_PUT_VAR(ncid,time_id    ,grid%XTIME)
    status=NF90_PUT_VAR(ncid,XTIME_id   ,grid%XTIME)
    status=NF90_PUT_VAR(ncid,XLAT_M_id  ,grid%XLAT_M    ,start=(/1,1,1/),count=(/grid%x_grd_num     ,grid%y_grd_num     ,1/))
    status=NF90_PUT_VAR(ncid,XLONG_M_id ,grid%XLONG_M   ,start=(/1,1,1/),count=(/grid%x_grd_num     ,grid%y_grd_num     ,1/))
    status=NF90_PUT_VAR(ncid,XLAT_U_id  ,grid%XLAT_U    ,start=(/1,1,1/),count=(/grid%x_grd_num_u   ,grid%y_grd_num_u   ,1/))
    status=NF90_PUT_VAR(ncid,XLONG_U_id ,grid%XLONG_U   ,start=(/1,1,1/),count=(/grid%x_grd_num_u   ,grid%y_grd_num_u   ,1/))
    status=NF90_PUT_VAR(ncid,XLAT_V_id  ,grid%XLAT_V    ,start=(/1,1,1/),count=(/grid%x_grd_num_v   ,grid%y_grd_num_v   ,1/))
    status=NF90_PUT_VAR(ncid,XLONG_V_id ,grid%XLONG_V   ,start=(/1,1,1/),count=(/grid%x_grd_num_v   ,grid%y_grd_num_v   ,1/))
    status=NF90_PUT_VAR(ncid,Z_id       ,grid%za        ,start=(/1,1,1/),count=(/grid%x_grd_num     ,grid%y_grd_num     ,1/))
    status=NF90_PUT_VAR(ncid,U_id       ,grid%ua        ,start=(/1,1,1/),count=(/grid%x_grd_num_u   ,grid%y_grd_num_u   ,1/))
    status=NF90_PUT_VAR(ncid,V_id       ,grid%va        ,start=(/1,1,1/),count=(/grid%x_grd_num_v   ,grid%y_grd_num_v   ,1/))

    status=NF90_CLOSE(ncid)
    end subroutine create_output



    !Write the output data into the netCDF file
    subroutine write_output(step,z,u,v)
    implicit none
    integer,intent(in)  :: step
    real,intent(in)     :: Z(grid%x_grd_num     ,grid%y_grd_num)    ,&
                           U(grid%x_grd_num_u   ,grid%y_grd_num_u)  ,&
                           V(grid%x_grd_num_v   ,grid%y_grd_num_v)

    !netCDF
    integer ncid,status

    !Variables ID
    integer XLAT_M_id  ,&
            XLONG_M_id ,&
            XLAT_U_id  ,&
            XLONG_U_id ,&
            XLAT_V_id  ,&
            XLONG_V_id ,&
            Z_id,&
            U_id,&
            V_id

    status=NF90_OPEN(trim(adjustl(output_data_path))//slash//'BPEM_output.nc',OR(NF90_WRITE,NF90_SHARE),ncid)
    status=NF90_INQ_VARID(ncid,'XLAT_M' ,XLAT_M_id  )
    status=NF90_INQ_VARID(ncid,'XLONG_M',XLONG_M_id )
    status=NF90_INQ_VARID(ncid,'XLAT_U' ,XLAT_U_id  )
    status=NF90_INQ_VARID(ncid,'XLONG_U',XLONG_U_id )
    status=NF90_INQ_VARID(ncid,'XLAT_V' ,XLAT_V_id  )
    status=NF90_INQ_VARID(ncid,'XLONG_V',XLONG_V_id )
    status=NF90_INQ_VARID(ncid,'Z',Z_id)
    status=NF90_INQ_VARID(ncid,'U',U_id)
    status=NF90_INQ_VARID(ncid,'V',V_id)
    status=NF90_PUT_VAR(ncid,XLAT_M_id  ,grid%XLAT_M    ,start=(/1,1,step/) ,count=(/grid%x_grd_num     ,grid%y_grd_num     ,1/))
    status=NF90_PUT_VAR(ncid,XLONG_M_id ,grid%XLONG_M   ,start=(/1,1,step/) ,count=(/grid%x_grd_num     ,grid%y_grd_num     ,1/))
    status=NF90_PUT_VAR(ncid,XLAT_U_id  ,grid%XLAT_U    ,start=(/1,1,step/) ,count=(/grid%x_grd_num_u   ,grid%y_grd_num_u   ,1/))
    status=NF90_PUT_VAR(ncid,XLONG_U_id ,grid%XLONG_U   ,start=(/1,1,step/) ,count=(/grid%x_grd_num_u   ,grid%y_grd_num_u   ,1/))
    status=NF90_PUT_VAR(ncid,XLAT_V_id  ,grid%XLAT_V    ,start=(/1,1,step/) ,count=(/grid%x_grd_num_v   ,grid%y_grd_num_v   ,1/))
    status=NF90_PUT_VAR(ncid,XLONG_V_id ,grid%XLONG_V   ,start=(/1,1,step/) ,count=(/grid%x_grd_num_v   ,grid%y_grd_num_v   ,1/))
    status=NF90_PUT_VAR(ncid,Z_id,Z                     ,start=(/1,1,step/) ,count=(/grid%x_grd_num     ,grid%y_grd_num     ,1/))
    status=NF90_PUT_VAR(ncid,U_id,U                     ,start=(/1,1,step/) ,count=(/grid%x_grd_num_u   ,grid%y_grd_num_u   ,1/))
    status=NF90_PUT_VAR(ncid,V_id,V                     ,start=(/1,1,step/) ,count=(/grid%x_grd_num_v   ,grid%y_grd_num_v   ,1/))
    status=NF90_CLOSE(ncid)

    end subroutine write_output

end module output