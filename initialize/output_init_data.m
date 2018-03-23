function output_init_data(output_path,grid)

mode=  netcdf.getConstant('NETCDF4');
mode = bitor(mode,netcdf.getConstant('CLOBBER'));
ncid=netcdf.create(output_path,mode);
compress_level=5;
disp(['ncid=',num2str(ncid)])

disp('## Defining Global Attributes...')
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'BDY_WIDTH'                  ,grid.bdy_width);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'BDY_TIME_NUM'               ,size(grid.bdy_times,2));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'BDY_TEND_TIME_NUM'          ,size(grid.bdy_times,2)-1);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'CEN_LAT'                    ,grid.CEN_LAT);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'WEST-EAST_GRID_DIMENSION'   ,grid.E_WE);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'SOUTH-NORTH_GRID_DIMENSION' ,grid.E_SN);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'STAND_LON'                  ,grid.STAND_LON);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'TRUELAT1'                   ,grid.TRUELAT1);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'TRUELAT2'                   ,grid.TRUELAT2);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'MAP_PROJ'                   ,grid.MAP_PROJ);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'CEN_LON'                    ,grid.CEN_LON);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'DX'                         ,grid.DX);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'DY'                         ,grid.DY);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'start_year'                 ,grid.start_year)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'start_month'                ,grid.start_month)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'start_day'                  ,grid.start_day)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'start_hour'                 ,grid.start_hour)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'start_minute'               ,grid.start_minute)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'start_second'               ,grid.start_second)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'end_year'                   ,grid.end_year)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'end_month'                  ,grid.end_month)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'end_day'                    ,grid.end_day)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'end_hour'                   ,grid.end_hour)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'end_minute'                 ,grid.end_minute)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'end_second'                 ,grid.end_second)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'interval_seconds'           ,grid.interval_seconds)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'spec_zone'                  ,grid.spec_zone)
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'relax_zone'                 ,grid.relax_zone)

west_east_dimID         = netcdf.defDim(ncid,'west_east',grid.E_WE-1);
south_north_dimID       = netcdf.defDim(ncid,'south_north',grid.E_SN-1);
west_east_stag_dimID    = netcdf.defDim(ncid,'west_east_stag',grid.E_WE);
south_north_stag_dimID  = netcdf.defDim(ncid,'south_north_stag',grid.E_SN);
bdy_time_dimID          = netcdf.defDim(ncid,'bdy_time',size(grid.bdy_times,2));
bdy_tend_time_dimID     = netcdf.defDim(ncid,'bdy_tend_time',size(grid.bdy_times,2)-1);
DateStrLen_dimID        = netcdf.defDim(ncid,'DateStrLen',size(grid.bdy_times,1));
bdy_width_dimID         = netcdf.defDim(ncid,'bdy_width',grid.bdy_width);

disp('## Defining Variables...')
bdy_times_ID = netcdf.defVar(ncid,'bdy_times','NC_CHAR',[DateStrLen_dimID,bdy_time_dimID]);
netcdf.defVarDeflate(ncid,bdy_times_ID,true,true,compress_level);

XLAT_M_id = netcdf.defVar(ncid,'XLAT_M','float',[west_east_dimID,south_north_dimID]);
netcdf.defVarDeflate(ncid,XLAT_M_id,true,true,compress_level);
netcdf.putAtt(ncid,XLAT_M_id,'units','degrees latitude');
netcdf.putAtt(ncid,XLAT_M_id,'description','Latitude on mass grid');

XLONG_M_id = netcdf.defVar(ncid,'XLONG_M','float',[west_east_dimID,south_north_dimID]);
netcdf.defVarDeflate(ncid,XLONG_M_id,true,true,compress_level);
netcdf.putAtt(ncid,XLONG_M_id,'units','degrees longitude');
netcdf.putAtt(ncid,XLONG_M_id,'description','Longitude on mass grid');

XLAT_C_id = netcdf.defVar(ncid,'XLAT_C','float',[west_east_stag_dimID,south_north_stag_dimID]);
netcdf.defVarDeflate(ncid,XLAT_C_id,true,true,compress_level);
netcdf.putAtt(ncid,XLAT_C_id,'units','degrees latitude');
netcdf.putAtt(ncid,XLAT_C_id,'description','Latitude at grid cell corners');

XLONG_C_id = netcdf.defVar(ncid,'XLONG_C','float',[west_east_stag_dimID,south_north_stag_dimID]);
netcdf.defVarDeflate(ncid,XLONG_C_id,true,true,compress_level);
netcdf.putAtt(ncid,XLONG_C_id,'units','degrees longitude');
netcdf.putAtt(ncid,XLONG_C_id,'description','Longitude at grid cell corners');

XLAT_V_id = netcdf.defVar(ncid,'XLAT_V','float',[west_east_dimID,south_north_stag_dimID]);
netcdf.defVarDeflate(ncid,XLAT_V_id,true,true,compress_level);
netcdf.putAtt(ncid,XLAT_V_id,'units','degrees latitude');
netcdf.putAtt(ncid,XLAT_V_id,'description','Latitude on V grid');

XLONG_V_id = netcdf.defVar(ncid,'XLONG_V','float',[west_east_dimID,south_north_stag_dimID]);
netcdf.defVarDeflate(ncid,XLONG_V_id,true,true,compress_level);
netcdf.putAtt(ncid,XLONG_V_id,'units','degrees longitude');
netcdf.putAtt(ncid,XLONG_V_id,'description','Longitude on V grid');

XLAT_U_id = netcdf.defVar(ncid,'XLAT_U','float',[west_east_stag_dimID,south_north_dimID]);
netcdf.defVarDeflate(ncid,XLAT_U_id,true,true,compress_level);
netcdf.putAtt(ncid,XLAT_U_id,'units','degrees latitude');
netcdf.putAtt(ncid,XLAT_U_id,'description','Latitude on U grid');

XLONG_U_id = netcdf.defVar(ncid,'XLONG_U','float',[west_east_stag_dimID,south_north_dimID]);
netcdf.defVarDeflate(ncid,XLONG_U_id,true,true,compress_level);
netcdf.putAtt(ncid,XLONG_U_id,'units','degrees longitude');
netcdf.putAtt(ncid,XLONG_U_id,'description','Longitude on U grid');

CLAT_id = netcdf.defVar(ncid,'CLAT','float',[west_east_dimID,south_north_dimID]);
netcdf.defVarDeflate(ncid,CLAT_id,true,true,compress_level);
netcdf.putAtt(ncid,CLAT_id,'units','degrees latitude');
netcdf.putAtt(ncid,CLAT_id,'description','Computational latitude on mass grid');

CLONG_id = netcdf.defVar(ncid,'CLONG','float',[west_east_dimID,south_north_dimID]);
netcdf.defVarDeflate(ncid,CLONG_id,true,true,compress_level);
netcdf.putAtt(ncid,CLONG_id,'units','degrees longitude');
netcdf.putAtt(ncid,CLONG_id,'description','Computational longitude on mass grid');

MAPFAC_M_id = netcdf.defVar(ncid,'MAPFAC_M','float',[west_east_dimID,south_north_dimID]);
netcdf.defVarDeflate(ncid,MAPFAC_M_id,true,true,compress_level);
netcdf.putAtt(ncid,MAPFAC_M_id,'units','none');
netcdf.putAtt(ncid,MAPFAC_M_id,'description','Mapfactor on mass grid');

MAPFAC_C_id = netcdf.defVar(ncid,'MAPFAC_C','float',[west_east_stag_dimID,south_north_stag_dimID]);
netcdf.defVarDeflate(ncid,MAPFAC_C_id,true,true,compress_level);
netcdf.putAtt(ncid,MAPFAC_C_id,'units','none');
netcdf.putAtt(ncid,MAPFAC_C_id,'description','Mapfactor on the cell conrners');

MAPFAC_V_id = netcdf.defVar(ncid,'MAPFAC_V','float',[west_east_dimID,south_north_stag_dimID]);
netcdf.defVarDeflate(ncid,MAPFAC_V_id,true,true,compress_level);
netcdf.putAtt(ncid,MAPFAC_V_id,'units','none');
netcdf.putAtt(ncid,MAPFAC_V_id,'description','Mapfactor on V grid');

MAPFAC_U_id = netcdf.defVar(ncid,'MAPFAC_U','float',[west_east_stag_dimID,south_north_dimID]);
netcdf.defVarDeflate(ncid,MAPFAC_U_id,true,true,compress_level);
netcdf.putAtt(ncid,MAPFAC_U_id,'units','none');
netcdf.putAtt(ncid,MAPFAC_U_id,'description','Mapfactor on U grid');

MAPFAC_MX_id = netcdf.defVar(ncid,'MAPFAC_MX','float',[west_east_dimID,south_north_dimID]);
netcdf.defVarDeflate(ncid,MAPFAC_MX_id,true,true,compress_level);
netcdf.putAtt(ncid,MAPFAC_MX_id,'units','none');
netcdf.putAtt(ncid,MAPFAC_MX_id,'description','Mapfactor (x-dir) on mass grid');

MAPFAC_VX_id = netcdf.defVar(ncid,'MAPFAC_VX','float',[west_east_dimID,south_north_stag_dimID]);
netcdf.defVarDeflate(ncid,MAPFAC_VX_id,true,true,compress_level);
netcdf.putAtt(ncid,MAPFAC_VX_id,'units','none');
netcdf.putAtt(ncid,MAPFAC_VX_id,'description','Mapfactor (x-dir) on V grid');

MAPFAC_UX_id = netcdf.defVar(ncid,'MAPFAC_UX','float',[west_east_stag_dimID,south_north_dimID]);
netcdf.defVarDeflate(ncid,MAPFAC_UX_id,true,true,compress_level);
netcdf.putAtt(ncid,MAPFAC_UX_id,'units','none');
netcdf.putAtt(ncid,MAPFAC_UX_id,'description','Mapfactor (x-dir) on U grid');

MAPFAC_MY_id = netcdf.defVar(ncid,'MAPFAC_MY','float',[west_east_dimID,south_north_dimID]);
netcdf.defVarDeflate(ncid,MAPFAC_MY_id,true,true,compress_level);
netcdf.putAtt(ncid,MAPFAC_MY_id,'units','none');
netcdf.putAtt(ncid,MAPFAC_MY_id,'description','Mapfactor (y-dir) on mass grid');

MAPFAC_VY_id = netcdf.defVar(ncid,'MAPFAC_VY','float',[west_east_dimID,south_north_stag_dimID]);
netcdf.defVarDeflate(ncid,MAPFAC_VY_id,true,true,compress_level);
netcdf.putAtt(ncid,MAPFAC_VY_id,'units','none');
netcdf.putAtt(ncid,MAPFAC_VY_id,'description','Mapfactor (y-dir) on V grid');

MAPFAC_UY_id = netcdf.defVar(ncid,'MAPFAC_UY','float',[west_east_stag_dimID,south_north_dimID]);
netcdf.defVarDeflate(ncid,MAPFAC_UY_id,true,true,compress_level);
netcdf.putAtt(ncid,MAPFAC_UY_id,'units','none');
netcdf.putAtt(ncid,MAPFAC_UY_id,'description','Mapfactor (y-dir) on U grid');

E_id = netcdf.defVar(ncid,'E','float',[west_east_dimID,south_north_dimID]);
netcdf.defVarDeflate(ncid,E_id,true,true,compress_level);
netcdf.putAtt(ncid,E_id,'units','-');
netcdf.putAtt(ncid,E_id,'description','Coriolis E parameter');

F_M_id = netcdf.defVar(ncid,'F_M','float',[west_east_dimID,south_north_dimID]);
netcdf.defVarDeflate(ncid,F_M_id,true,true,compress_level);
netcdf.putAtt(ncid,F_M_id,'units','s^-1');
netcdf.putAtt(ncid,F_M_id,'description','Coriolis F parameter for mass grid');

F_C_id = netcdf.defVar(ncid,'F_C','float',[west_east_stag_dimID,south_north_stag_dimID]);
netcdf.defVarDeflate(ncid,F_C_id,true,true,compress_level);
netcdf.putAtt(ncid,F_C_id,'units','s^-1');
netcdf.putAtt(ncid,F_C_id,'description','Coriolis F parameter for the cell conrners');

F_U_id = netcdf.defVar(ncid,'F_U','float',[west_east_stag_dimID,south_north_dimID]);
netcdf.defVarDeflate(ncid,F_U_id,true,true,compress_level);
netcdf.putAtt(ncid,F_U_id,'units','s^-1');
netcdf.putAtt(ncid,F_U_id,'description','Coriolis F parameter for U grid');

F_V_id = netcdf.defVar(ncid,'F_V','float',[west_east_dimID,south_north_stag_dimID]);
netcdf.defVarDeflate(ncid,F_V_id,true,true,compress_level);
netcdf.putAtt(ncid,F_V_id,'units','s^-1');
netcdf.putAtt(ncid,F_V_id,'description','Coriolis F parameter for V grid');

Z_id = netcdf.defVar(ncid,'Z','float',[west_east_dimID,south_north_dimID,bdy_time_dimID]);
netcdf.defVarDeflate(ncid,Z_id,true,true,compress_level);
netcdf.putAtt(ncid,Z_id,'units','dagpm');
netcdf.putAtt(ncid,Z_id,'description','Geopotential Height');

U_id = netcdf.defVar(ncid,'U','float',[west_east_stag_dimID,south_north_dimID,bdy_time_dimID]);
netcdf.defVarDeflate(ncid,U_id,true,true,compress_level);
netcdf.putAtt(ncid,U_id,'units','m/s');
netcdf.putAtt(ncid,U_id,'description','U-Component Wind');

V_id = netcdf.defVar(ncid,'V','float',[west_east_dimID,south_north_stag_dimID,bdy_time_dimID]);
netcdf.defVarDeflate(ncid,V_id,true,true,compress_level);
netcdf.putAtt(ncid,V_id,'units','m/s');
netcdf.putAtt(ncid,V_id,'description','V-Component Wind');

Z_BYS_id = netcdf.defVar(ncid,'Z_BYS','float',[west_east_dimID,bdy_width_dimID,bdy_time_dimID]);
netcdf.defVarDeflate(ncid,Z_BYS_id,true,true,compress_level);
netcdf.putAtt(ncid,Z_BYS_id,'units','dagpm');
netcdf.putAtt(ncid,Z_BYS_id,'description','Geopotential Height bdy at Y-axis start');

Z_BYE_id = netcdf.defVar(ncid,'Z_BYE','float',[west_east_dimID,bdy_width_dimID,bdy_time_dimID]);
netcdf.defVarDeflate(ncid,Z_BYE_id,true,true,compress_level);
netcdf.putAtt(ncid,Z_BYE_id,'units','dagpm');
netcdf.putAtt(ncid,Z_BYE_id,'description','Geopotential Height bdy at Y-axis end');

Z_BXS_id = netcdf.defVar(ncid,'Z_BXS','float',[south_north_dimID,bdy_width_dimID,bdy_time_dimID]);
netcdf.defVarDeflate(ncid,Z_BXS_id,true,true,compress_level);
netcdf.putAtt(ncid,Z_BXS_id,'units','dagpm');
netcdf.putAtt(ncid,Z_BXS_id,'description','Geopotential Height bdy at X-axis start');

Z_BXE_id = netcdf.defVar(ncid,'Z_BXE','float',[south_north_dimID,bdy_width_dimID,bdy_time_dimID]);
netcdf.defVarDeflate(ncid,Z_BXE_id,true,true,compress_level);
netcdf.putAtt(ncid,Z_BXE_id,'units','dagpm');
netcdf.putAtt(ncid,Z_BXE_id,'description','Geopotential Height bdy at X-axis end');

U_BYS_id = netcdf.defVar(ncid,'U_BYS','float',[west_east_stag_dimID,bdy_width_dimID,bdy_time_dimID]);
netcdf.defVarDeflate(ncid,U_BYS_id,true,true,compress_level);
netcdf.putAtt(ncid,U_BYS_id,'units','m/s');
netcdf.putAtt(ncid,U_BYS_id,'description','U-Component Wind bdy at Y-axis start');

U_BYE_id = netcdf.defVar(ncid,'U_BYE','float',[west_east_stag_dimID,bdy_width_dimID,bdy_time_dimID]);
netcdf.defVarDeflate(ncid,U_BYE_id,true,true,compress_level);
netcdf.putAtt(ncid,U_BYE_id,'units','m/s');
netcdf.putAtt(ncid,U_BYE_id,'description','U-Component Wind bdy at Y-axis end');

U_BXS_id = netcdf.defVar(ncid,'U_BXS','float',[south_north_dimID,bdy_width_dimID,bdy_time_dimID]);
netcdf.defVarDeflate(ncid,U_BXS_id,true,true,compress_level);
netcdf.putAtt(ncid,U_BXS_id,'units','m/s');
netcdf.putAtt(ncid,U_BXS_id,'description','U-Component Wind bdy at X-axis start');

U_BXE_id = netcdf.defVar(ncid,'U_BXE','float',[south_north_dimID,bdy_width_dimID,bdy_time_dimID]);
netcdf.defVarDeflate(ncid,U_BXE_id,true,true,compress_level);
netcdf.putAtt(ncid,U_BXE_id,'units','m/s');
netcdf.putAtt(ncid,U_BXE_id,'description','U-Component Wind bdy at X-axis end');

V_BYS_id = netcdf.defVar(ncid,'V_BYS','float',[west_east_dimID,bdy_width_dimID,bdy_time_dimID]);
netcdf.defVarDeflate(ncid,V_BYS_id,true,true,compress_level);
netcdf.putAtt(ncid,V_BYS_id,'units','m/s');
netcdf.putAtt(ncid,V_BYS_id,'description','V-Component Wind bdy at Y-axis start');

V_BYE_id = netcdf.defVar(ncid,'V_BYE','float',[west_east_dimID,bdy_width_dimID,bdy_time_dimID]);
netcdf.defVarDeflate(ncid,V_BYE_id,true,true,compress_level);
netcdf.putAtt(ncid,V_BYE_id,'units','m/s');
netcdf.putAtt(ncid,V_BYE_id,'description','V-Component Wind bdy at Y-axis end');

V_BXS_id = netcdf.defVar(ncid,'V_BXS','float',[south_north_stag_dimID,bdy_width_dimID,bdy_time_dimID]);
netcdf.defVarDeflate(ncid,V_BXS_id,true,true,compress_level);
netcdf.putAtt(ncid,V_BXS_id,'units','m/s');
netcdf.putAtt(ncid,V_BXS_id,'description','V-Component Wind bdy at X-axis start');

V_BXE_id = netcdf.defVar(ncid,'V_BXE','float',[south_north_stag_dimID,bdy_width_dimID,bdy_time_dimID]);
netcdf.defVarDeflate(ncid,V_BXE_id,true,true,compress_level);
netcdf.putAtt(ncid,V_BXE_id,'units','m/s');
netcdf.putAtt(ncid,V_BXE_id,'description','V-Component Wind bdy at X-axis end');

Z_BTYS_id = netcdf.defVar(ncid,'Z_BTYS','float',[west_east_dimID,bdy_width_dimID,bdy_tend_time_dimID]);
netcdf.defVarDeflate(ncid,Z_BTYS_id,true,true,compress_level);
netcdf.putAtt(ncid,Z_BTYS_id,'units','dagpm');
netcdf.putAtt(ncid,Z_BTYS_id,'description','Geopotential Height bdy tend at Y-axis start');

Z_BTYE_id = netcdf.defVar(ncid,'Z_BTYE','float',[west_east_dimID,bdy_width_dimID,bdy_tend_time_dimID]);
netcdf.defVarDeflate(ncid,Z_BTYE_id,true,true,compress_level);
netcdf.putAtt(ncid,Z_BTYE_id,'units','dagpm');
netcdf.putAtt(ncid,Z_BTYE_id,'description','Geopotential Height bdy tend at Y-axis end');

Z_BTXS_id = netcdf.defVar(ncid,'Z_BTXS','float',[south_north_dimID,bdy_width_dimID,bdy_tend_time_dimID]);
netcdf.defVarDeflate(ncid,Z_BTXS_id,true,true,compress_level);
netcdf.putAtt(ncid,Z_BTXS_id,'units','dagpm');
netcdf.putAtt(ncid,Z_BTXS_id,'description','Geopotential Height bdy tend at X-axis start');

Z_BTXE_id = netcdf.defVar(ncid,'Z_BTXE','float',[south_north_dimID,bdy_width_dimID,bdy_tend_time_dimID]);
netcdf.defVarDeflate(ncid,Z_BTXE_id,true,true,compress_level);
netcdf.putAtt(ncid,Z_BTXE_id,'units','dagpm');
netcdf.putAtt(ncid,Z_BTXE_id,'description','Geopotential Height bdy tend at X-axis end');

U_BTYS_id = netcdf.defVar(ncid,'U_BTYS','float',[west_east_stag_dimID,bdy_width_dimID,bdy_tend_time_dimID]);
netcdf.defVarDeflate(ncid,U_BTYS_id,true,true,compress_level);
netcdf.putAtt(ncid,U_BTYS_id,'units','m/s');
netcdf.putAtt(ncid,U_BTYS_id,'description','U-Component Wind bdy tend at Y-axis start');

U_BTYE_id = netcdf.defVar(ncid,'U_BTYE','float',[west_east_stag_dimID,bdy_width_dimID,bdy_tend_time_dimID]);
netcdf.defVarDeflate(ncid,U_BTYE_id,true,true,compress_level);
netcdf.putAtt(ncid,U_BTYE_id,'units','m/s');
netcdf.putAtt(ncid,U_BTYE_id,'description','U-Component Wind bdy tend at Y-axis end');

U_BTXS_id = netcdf.defVar(ncid,'U_BTXS','float',[south_north_dimID,bdy_width_dimID,bdy_tend_time_dimID]);
netcdf.defVarDeflate(ncid,U_BTXS_id,true,true,compress_level);
netcdf.putAtt(ncid,U_BTXS_id,'units','m/s');
netcdf.putAtt(ncid,U_BTXS_id,'description','U-Component Wind bdy tend at X-axis start');

U_BTXE_id = netcdf.defVar(ncid,'U_BTXE','float',[south_north_dimID,bdy_width_dimID,bdy_tend_time_dimID]);
netcdf.defVarDeflate(ncid,U_BTXE_id,true,true,compress_level);
netcdf.putAtt(ncid,U_BTXE_id,'units','m/s');
netcdf.putAtt(ncid,U_BTXE_id,'description','U-Component Wind bdy tend at X-axis end');

V_BTYS_id = netcdf.defVar(ncid,'V_BTYS','float',[west_east_dimID,bdy_width_dimID,bdy_tend_time_dimID]);
netcdf.defVarDeflate(ncid,V_BTYS_id,true,true,compress_level);
netcdf.putAtt(ncid,V_BTYS_id,'units','m/s');
netcdf.putAtt(ncid,V_BTYS_id,'description','V-Component Wind bdy tend at Y-axis start');

V_BTYE_id = netcdf.defVar(ncid,'V_BTYE','float',[west_east_dimID,bdy_width_dimID,bdy_tend_time_dimID]);
netcdf.defVarDeflate(ncid,V_BTYE_id,true,true,compress_level);
netcdf.putAtt(ncid,V_BTYE_id,'units','m/s');
netcdf.putAtt(ncid,V_BTYE_id,'description','V-Component Wind bdy tend at Y-axis end');

V_BTXS_id = netcdf.defVar(ncid,'V_BTXS','float',[south_north_stag_dimID,bdy_width_dimID,bdy_tend_time_dimID]);
netcdf.defVarDeflate(ncid,V_BTXS_id,true,true,compress_level);
netcdf.putAtt(ncid,V_BTXS_id,'units','m/s');
netcdf.putAtt(ncid,V_BTXS_id,'description','V-Component Wind bdy tend at X-axis start');

V_BTXE_id = netcdf.defVar(ncid,'V_BTXE','float',[south_north_stag_dimID,bdy_width_dimID,bdy_tend_time_dimID]);
netcdf.defVarDeflate(ncid,V_BTXE_id,true,true,compress_level);
netcdf.putAtt(ncid,V_BTXE_id,'units','m/s');
netcdf.putAtt(ncid,V_BTXE_id,'description','V-Component Wind bdy tend at X-axis end');

disp('## Putting Variables...')
netcdf.putVar(ncid, bdy_times_ID,grid.bdy_times);
netcdf.putVar(ncid, XLONG_M_id,  grid.XLONG_M);
netcdf.putVar(ncid, XLAT_M_id,   grid.XLAT_M);
netcdf.putVar(ncid, XLONG_U_id,  grid.XLONG_U);
netcdf.putVar(ncid, XLAT_U_id,   grid.XLAT_U);
netcdf.putVar(ncid, XLONG_V_id,  grid.XLONG_V);
netcdf.putVar(ncid, XLAT_V_id,   grid.XLAT_V);
netcdf.putVar(ncid, XLONG_C_id,  grid.XLONG_C);
netcdf.putVar(ncid, XLAT_C_id,   grid.XLAT_C);
netcdf.putVar(ncid, CLONG_id,    grid.CLONG);
netcdf.putVar(ncid, CLAT_id,     grid.CLAT);
netcdf.putVar(ncid, MAPFAC_M_id, grid.MAPFAC_M);
netcdf.putVar(ncid, MAPFAC_C_id, grid.MAPFAC_C);
netcdf.putVar(ncid, MAPFAC_V_id, grid.MAPFAC_V);
netcdf.putVar(ncid, MAPFAC_U_id, grid.MAPFAC_U);
netcdf.putVar(ncid, MAPFAC_MX_id,grid.MAPFAC_MX);
netcdf.putVar(ncid, MAPFAC_VX_id,grid.MAPFAC_VX);
netcdf.putVar(ncid, MAPFAC_UX_id,grid.MAPFAC_UX);
netcdf.putVar(ncid, MAPFAC_MY_id,grid.MAPFAC_MY);
netcdf.putVar(ncid, MAPFAC_VY_id,grid.MAPFAC_VY);
netcdf.putVar(ncid, MAPFAC_UY_id,grid.MAPFAC_UY);
netcdf.putVar(ncid, E_id,        grid.E);
netcdf.putVar(ncid, F_M_id,      grid.F_M);
netcdf.putVar(ncid, F_C_id,      grid.F_C);
netcdf.putVar(ncid, F_U_id,      grid.F_U);
netcdf.putVar(ncid, F_V_id,      grid.F_V);
netcdf.putVar(ncid, Z_id,        grid.Z);
netcdf.putVar(ncid, U_id,        grid.U);
netcdf.putVar(ncid, V_id,        grid.V);

netcdf.putVar(ncid, Z_BYS_id,    grid.Z_BYS);
netcdf.putVar(ncid, Z_BYE_id,    grid.Z_BYE);
netcdf.putVar(ncid, Z_BXS_id,    grid.Z_BXS);
netcdf.putVar(ncid, Z_BXE_id,    grid.Z_BXE);

netcdf.putVar(ncid, U_BYS_id,    grid.U_BYS);
netcdf.putVar(ncid, U_BYE_id,    grid.U_BYE);
netcdf.putVar(ncid, U_BXS_id,    grid.U_BXS);
netcdf.putVar(ncid, U_BXE_id,    grid.U_BXE);

netcdf.putVar(ncid, V_BYS_id,    grid.V_BYS);
netcdf.putVar(ncid, V_BYE_id,    grid.V_BYE);
netcdf.putVar(ncid, V_BXS_id,    grid.V_BXS);
netcdf.putVar(ncid, V_BXE_id,    grid.V_BXE);

netcdf.putVar(ncid, Z_BTYS_id,   grid.Z_BTYS);
netcdf.putVar(ncid, Z_BTYE_id,   grid.Z_BTYE);
netcdf.putVar(ncid, Z_BTXS_id,   grid.Z_BTXS);
netcdf.putVar(ncid, Z_BTXE_id,   grid.Z_BTXE);

netcdf.putVar(ncid, U_BTYS_id,   grid.U_BTYS);
netcdf.putVar(ncid, U_BTYE_id,   grid.U_BTYE);
netcdf.putVar(ncid, U_BTXS_id,   grid.U_BTXS);
netcdf.putVar(ncid, U_BTXE_id,	 grid.U_BTXE);

netcdf.putVar(ncid, V_BTYS_id,   grid.V_BTYS);
netcdf.putVar(ncid, V_BTYE_id,   grid.V_BTYE);
netcdf.putVar(ncid, V_BTXS_id,   grid.V_BTXS);
netcdf.putVar(ncid, V_BTXE_id,   grid.V_BTXE);

netcdf.close(ncid)