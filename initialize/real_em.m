% Written by Zhou Lilong at Dec 2017
% This program is used to create the initial field for BPEM
% 2 input files are required by this program
% 
% source_nc : a ERA-interim netCDF file, which should include "longitude",
% "latitude", "z","u","v" filed. For an original atmosphere barotropic
% model, the data on 500hPa is original choice.
% 
% grid_nc : a geo_em file, which created by WRF Preprocess System(WPS)
% 
% use_WRF_geo_em : when set to 1, this program will try to obtain the geo
% info form WRF geo_em.d01.nc file. Otherwise, you need to set the
% parameters for the lambert projection, and the domain info will be
% generate by your input.
% 
% input_path : where to find the source_nc file

% clc
% clear
t1=clock;
%%%%%%%%%%%%
% Settings %
%%%%%%%%%%%%
% all of the variables shall be set as integer %
start_year      = 2017;
start_month     = 10;
start_day       = 01;
start_hour      = 00;
start_minute    = 00;
start_second    = 00;
end_year        = 2017;
end_month       = 10;
end_day         = 31;
end_hour        = 00;
end_minute      = 00;
end_second      = 00;
interval_seconds= 21600;
bdy_width       = 8;
spec_zone       = 2;
relax_zone      = bdy_width-spec_zone; %Do not modify!

input_path      = 'E:\study\BPEM\input';
source_nc       = '201710_ZUV.nc';

use_WRF_geo_em  = 0;
grid_nc         = 'geo_em.d01.nc';
%%%%%%%%%%%%%%%%%%%%%%%
% Parameters' setting %
%%%%%%%%%%%%%%%%%%%%%%%
% Lambert projection
%Input
%ref_lat	% the latitude of the center point 
%ref_lon	% the longitude of the center point
%dx         % the grid distance in the x-direction
%e_we       % the grid number since x direction
%e_sn       % the grid number since y direction
%stand_lon	% the longitude parallel with the y-axis
%truelat1	% the low true latitude for the Lambert 
%truelat2	% the high true latitude for the Lambert
ref_lat	 = 37.;
ref_lon	 = 106.;
dx       = 100000.;
e_we     = 100;
e_sn     = 80;
stand_lon= 106.;
truelat1 = 37.;
truelat2 = 37.;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set slash direction in the file system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
system='Windows'; % Choose from Windows or Linux
if strcmp(system,'Windows')
    slash='\';
elseif strcmp(system,'Linux')
    slash='/';
end

%%%%%%%%%%%%%%%%%%%
% Executable code %
%%%%%%%%%%%%%%%%%%%
input_source        = [input_path,slash,source_nc];
input_grid          = [input_path,slash,grid_nc];

source_time         = double(ncread(input_source,'time'));
interval_hours      = source_time(2)-source_time(1);
timenum0            = datenum('19000101000000','yyyymmddHHMMSS');
source_time         = source_time/24+timenum0;

source_time_temp    = str2num(datestr(source_time,'yyyymmddHHMMSS'));
start_time          = start_year*10^10+start_month*10^8+start_day*10^6+start_hour*10^4+start_minute*10^2+start_second;
end_time            = end_year*10^10+end_month*10^8+end_day*10^6+end_hour*10^4+end_minute*10^2+end_second;
start_time_index    = find(source_time_temp==start_time);
end_time_index      = find(source_time_temp==end_time);
bdy_time_num        = end_time_index-start_time_index+1;

source_longitude    = double(ncread(input_source,'longitude'));
source_latitude     = double(ncread(input_source,'latitude'));
source_longitude    = repmat(source_longitude,1,size(source_latitude,1));
source_latitude     = repmat(source_latitude,1,size(source_longitude,1));
source_latitude     = source_latitude';

source_z            = double(ncread(input_source,'z'));
source_u            = double(ncread(input_source,'u'));
source_v            = double(ncread(input_source,'v'));

%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the geo info %
%%%%%%%%%%%%%%%%%%%%%%%
if use_WRF_geo_em
    % Read info from geo_em.d01.nc
    grid                = read_geo_em_info(input_grid);
else
    grid                = generate_domain(ref_lat,ref_lon,dx,e_we,e_sn,stand_lon,truelat1,truelat2);
end
grid.bdy_times      = datestr(source_time(start_time_index:end_time_index),'yyyy-mm-dd_HH:MM:SS');
grid.bdy_times      = grid.bdy_times';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate the source data to model grid %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid_Z          = zeros(size(grid.XLONG_M,1),size(grid.XLONG_M,2),bdy_time_num);
grid_U          = zeros(size(grid.XLONG_U,1),size(grid.XLONG_U,2),bdy_time_num);
grid_V          = zeros(size(grid.XLONG_V,1),size(grid.XLONG_V,2),bdy_time_num);
grid_XLONG_M    = grid.XLONG_M;
grid_XLAT_M     = grid.XLAT_M;
grid_XLONG_U    = grid.XLONG_U;
grid_XLAT_U     = grid.XLAT_U;
grid_XLONG_V    = grid.XLONG_V;
grid_XLAT_V     = grid.XLAT_V;

grid_XLONG_M(grid_XLONG_M<0)=360+grid_XLONG_M(grid_XLONG_M<0);
grid_XLONG_U(grid_XLONG_U<0)=360+grid_XLONG_U(grid_XLONG_U<0);
grid_XLONG_V(grid_XLONG_V<0)=360+grid_XLONG_V(grid_XLONG_V<0);

parfor i=1:bdy_time_num
    grid_Z(:,:,i)	= griddata(source_longitude,source_latitude,source_z(:,:,start_time_index+i-1),grid_XLONG_M,grid_XLAT_M,'cubic');
    grid_U(:,:,i)	= griddata(source_longitude,source_latitude,source_u(:,:,start_time_index+i-1),grid_XLONG_U,grid_XLAT_U,'cubic');
	grid_V(:,:,i)	= griddata(source_longitude,source_latitude,source_v(:,:,start_time_index+i-1),grid_XLONG_V,grid_XLAT_V,'cubic');
end
grid.Z=grid_Z;
grid.U=grid_U;
grid.V=grid_V;
clear grid_z grid_u grid_v grid_XLONG_M grid_XLAT_M grid_XLONG_U grid_XLAT_U grid_XLONG_V grid_XLAT_V
% Convert z to Geopotential Height
grid.Z = grid.Z/9.8;
% Calculate coliolis f parameters for u,v and conrers
omega   = 0.00007292115;
grid.F_U=2.*omega*sind(grid.XLAT_U);
grid.F_V=2.*omega*sind(grid.XLAT_V);
grid.F_C=2.*omega*sind(grid.XLAT_C);
grid.F_M=grid.F;
% Calculate map factor on corner grid
theta1=90-grid.TRUELAT1;
theta2=90-grid.TRUELAT2;
thetac=90-grid.XLAT_C;
if grid.TRUELAT1-grid.TRUELAT2<0.01
    cone_paramter=sind(grid.TRUELAT1);
else
    cone_paramter=(ln(sind(theta1)-ln(theta2)))/(ln(tand((theta1)/2))-ln(tand((theta2)/2)));
end
grid.MAPFAC_C=sind(theta1)./sind(thetac).*power(tand(thetac./2)./tand(theta1./2),cone_paramter);

%%%%%%%%%%%%%%%%%%%%%
% Generate Boundary %
%%%%%%%%%%%%%%%%%%%%%
grid.bdy_width=bdy_width;
grid.Z_BYS=grid.Z(:,1:grid.bdy_width,:);
grid.Z_BYE=fliplr(grid.Z(:,end-grid.bdy_width+1:end,:));
grid.Z_BXS=grid.Z(1:grid.bdy_width,:,:);
grid.Z_BXE=flipud(grid.Z(end-grid.bdy_width+1:end,:,:));

grid.U_BYS=grid.U(:,1:grid.bdy_width,:);
grid.U_BYE=fliplr(grid.U(:,end-grid.bdy_width+1:end,:));
grid.U_BXS=grid.U(1:grid.bdy_width,:,:);
grid.U_BXE=flipud(grid.U(end-grid.bdy_width+1:end,:,:));

grid.V_BYS=grid.V(:,1:grid.bdy_width,:);
grid.V_BYE=fliplr(grid.V(:,end-grid.bdy_width+1:end,:));
grid.V_BXS=grid.V(1:grid.bdy_width,:,:);
grid.V_BXE=flipud(grid.V(end-grid.bdy_width+1:end,:,:));

grid.Z_BXS=permute(grid.Z_BXS,[2,1,3]);
grid.Z_BXE=permute(grid.Z_BXE,[2,1,3]);
grid.U_BXS=permute(grid.U_BXS,[2,1,3]);
grid.U_BXE=permute(grid.U_BXE,[2,1,3]);
grid.V_BXS=permute(grid.V_BXS,[2,1,3]);
grid.V_BXE=permute(grid.V_BXE,[2,1,3]);

grid.Z_BTXS=(grid.Z_BXS(:,:,2:end)-grid.Z_BXS(:,:,1:end-1))./interval_hours;
grid.Z_BTXE=(grid.Z_BXE(:,:,2:end)-grid.Z_BXE(:,:,1:end-1))./interval_hours;
grid.U_BTXS=(grid.U_BXS(:,:,2:end)-grid.U_BXS(:,:,1:end-1))./interval_hours;
grid.U_BTXE=(grid.U_BXE(:,:,2:end)-grid.U_BXE(:,:,1:end-1))./interval_hours;
grid.V_BTXS=(grid.V_BXS(:,:,2:end)-grid.V_BXS(:,:,1:end-1))./interval_hours;
grid.V_BTXE=(grid.V_BXE(:,:,2:end)-grid.V_BXE(:,:,1:end-1))./interval_hours;

grid.Z_BTYS=(grid.Z_BYS(:,:,2:end)-grid.Z_BYS(:,:,1:end-1))./interval_hours;
grid.Z_BTYE=(grid.Z_BYE(:,:,2:end)-grid.Z_BYE(:,:,1:end-1))./interval_hours;
grid.U_BTYS=(grid.U_BYS(:,:,2:end)-grid.U_BYS(:,:,1:end-1))./interval_hours;
grid.U_BTYE=(grid.U_BYE(:,:,2:end)-grid.U_BYE(:,:,1:end-1))./interval_hours;
grid.V_BTYS=(grid.V_BYS(:,:,2:end)-grid.V_BYS(:,:,1:end-1))./interval_hours;
grid.V_BTYE=(grid.V_BYE(:,:,2:end)-grid.V_BYE(:,:,1:end-1))./interval_hours;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the initial & boundary data into the netCDF files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid.start_year         = start_year;
grid.start_month        = start_month;
grid.start_day          = start_day;
grid.start_hour         = start_hour;
grid.start_minute       = start_minute;
grid.start_second       = start_second;
grid.end_year           = end_year;
grid.end_month          = end_month;
grid.end_day            = end_day;
grid.end_hour           = end_hour;
grid.end_minute         = end_minute;
grid.end_second         = end_second;
grid.interval_seconds   = interval_seconds;
grid.spec_zone          = spec_zone;
grid.relax_zone         = relax_zone;

output_init_data('BPEM_input.nc',grid)

%%%%%%%%
% Plot %
%%%%%%%%
ChinaL=shaperead('D:\MeteoInfo\map\bou2_4p.shp');
bou2_4px=[ChinaL(:).X];
bou2_4py=[ChinaL(:).Y];
figure
grid.XLONG_M(grid.XLONG_M<0)=360+grid.XLONG_M(grid.XLONG_M<0);
plot_X_interval   = max(grid.CEN_LON-min(min(grid.XLONG_M)),max(max(grid.XLONG_M))-grid.CEN_LON);
plot_Y_interval   = max(grid.CEN_LAT-min(min(grid.XLAT_M)),max(max(grid.XLAT_M))-grid.CEN_LAT);
plot_right_lon      = grid.CEN_LON+plot_X_interval;
plot_left_lon       = grid.CEN_LON-plot_X_interval;
plot_upper_lat      = grid.CEN_LAT+plot_Y_interval;
plot_bottom_lat     = grid.CEN_LAT-plot_Y_interval;
m_proj('lambert','lon',[plot_left_lon,plot_right_lon],'lat',[plot_bottom_lat,plot_upper_lat]);
m_pcolor(grid.XLONG_M,grid.XLAT_M,grid.Z(:,:,2));
hold on
m_plot(bou2_4px,bou2_4py,'k');
m_grid;
m_coast('color','k','linewidth',1);
colormap(jet)
shading interp
colorbar

title('Geopotential Height at 500hPa','fontsize',10,'fontweight','bold');
print('-dpng','-zbuffer','-r300','z.png');

figure
grid.XLONG_U(grid.XLONG_U<0)=360+grid.XLONG_U(grid.XLONG_U<0);
plot_X_interval   = max(grid.CEN_LON-min(min(grid.XLONG_U)),max(max(grid.XLONG_U))-grid.CEN_LON);
plot_Y_interval   = max(grid.CEN_LAT-min(min(grid.XLAT_U)),max(max(grid.XLAT_U))-grid.CEN_LAT);
plot_right_lon      = grid.CEN_LON+plot_X_interval;
plot_left_lon       = grid.CEN_LON-plot_X_interval;
plot_upper_lat      = grid.CEN_LAT+plot_Y_interval;
plot_bottom_lat     = grid.CEN_LAT-plot_Y_interval;
m_proj('lambert','lon',[plot_left_lon,plot_right_lon],'lat',[plot_bottom_lat,plot_upper_lat]);
m_pcolor(grid.XLONG_U,grid.XLAT_U,grid.U(:,:,2));
hold on
m_plot(bou2_4px,bou2_4py,'k');
m_grid;
m_coast('color','k','linewidth',1);
colormap(jet)
shading interp
colorbar

title('U-Wind Component at 500hPa','fontsize',10,'fontweight','bold');
print('-dpng','-zbuffer','-r300','u.png');

figure
grid.XLONG_V(grid.XLONG_V<0)=360+grid.XLONG_V(grid.XLONG_V<0);
plot_X_interval   = max(grid.CEN_LON-min(min(grid.XLONG_V)),max(max(grid.XLONG_V))-grid.CEN_LON);
plot_Y_interval   = max(grid.CEN_LAT-min(min(grid.XLAT_V)),max(max(grid.XLAT_V))-grid.CEN_LAT);
plot_right_lon      = grid.CEN_LON+plot_X_interval;
plot_left_lon       = grid.CEN_LON-plot_X_interval;
plot_upper_lat      = grid.CEN_LAT+plot_Y_interval;
plot_bottom_lat     = grid.CEN_LAT-plot_Y_interval;
m_proj('lambert','lon',[plot_left_lon,plot_right_lon],'lat',[plot_bottom_lat,plot_upper_lat]);
m_pcolor(grid.XLONG_V,grid.XLAT_V,grid.V(:,:,2))
hold on
m_plot(bou2_4px,bou2_4py,'k');
m_grid;
m_coast('color','k','linewidth',1);
colormap(jet)
shading interp
colorbar

title('V-Wind Component at 500hPa','fontsize',10,'fontweight','bold');
print('-dpng','-zbuffer','-r300','v.png');


t2=clock;
disp(['It took ',num2str(etime(t2,t1)),' seconds to run this program'])