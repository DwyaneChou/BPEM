% Written by Zhou Lilong at Mar 2018
% This program is used to generate the domain info for BPEM
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
function [grid]=generate_domain(ref_lat,ref_lon,dx,e_we,e_sn,stand_lon,...
                                truelat1,truelat2)
% Default parameters, do not modify %
s_sn              = 1;
s_we              = 1;
n                 = 1;
% parent_id         = 1;
i_parent_start    = 1;
j_parent_start    = 1;
parent_grid_ratio = 1;

[lon_m,lat_m,lon_u,lat_u,lon_v,lat_v,...
 mapfac_m,mapfac_u,mapfac_v]=lambert(n,parent_grid_ratio,i_parent_start,j_parent_start,...
                                     ref_lat,ref_lon,dx,s_we,e_we,s_sn,e_sn,stand_lon,...
                                     truelat1,truelat2);
grid.XLONG_M   = lon_m;
grid.XLAT_M    = lat_m;
grid.XLONG_U   = lon_u;
grid.XLAT_U    = lat_u;
grid.XLONG_V   = lon_v;
grid.XLAT_V    = lat_v;

grid.MAP_PROJ  = 1;
grid.CEN_LON   = ref_lon;
grid.CEN_LAT   = ref_lat;
grid.DX        = dx;
grid.DY        = dx;
grid.E_WE      = e_we;
grid.E_SN      = e_sn;
grid.STAND_LON = stand_lon;
grid.TRUELAT1  = truelat1;
grid.TRUELAT2  = truelat2;

grid.MAPFAC_M  = mapfac_m;
grid.MAPFAC_U  = mapfac_u;
grid.MAPFAC_V  = mapfac_v;

omega          = 0.00007292115;
grid.F         = 2.*omega*sind(grid.XLAT_M);


grid.XLONG_C   = zeros(e_we,e_sn);
grid.XLAT_C    = zeros(e_we,e_sn);
grid.E         = zeros(size(grid.XLONG_M));
grid.CLONG     = zeros(size(grid.XLONG_M));
grid.CLAT      = zeros(size(grid.XLONG_M));
grid.MAPFAC_MX = zeros(size(grid.XLONG_M));
grid.MAPFAC_VX = zeros(size(grid.XLONG_V));
grid.MAPFAC_UX = zeros(size(grid.XLONG_U));
grid.MAPFAC_MY = zeros(size(grid.XLONG_M));
grid.MAPFAC_VY = zeros(size(grid.XLONG_V));
grid.MAPFAC_UY = zeros(size(grid.XLONG_U));