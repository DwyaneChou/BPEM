%Input
%n          % domain number
%parent_id(n) %n elements array
%parent_grid_ratio(n) % n elements array
%i_parent_start(n) % n elements array
%j_parent_start(n) % n elements array
%ref_lat	% the latitude of the center point 
%ref_lon	% the longitude of the center point
%dx         % the grid distance in the x-direction
%s_we=1     %default value is 1
%e_we
%s_sn=1     %default value is 1
%e_sn
%stand_lon	% the longitude parallel with the y-axis
%truelat1	% the low true latitude for the Lambert 
%truelat2	% the high true latitude for the Lambert
function [lon_m,lat_m,lon_u,lat_u,lon_v,lat_v,...
          mapfac_m,mapfac_u,mapfac_v]=lambert(n,parent_grid_ratio,i_parent_start,j_parent_start,...
                                              ref_lat,ref_lon,dx,s_we,e_we,s_sn,e_sn,stand_lon,truelat1,truelat2)


%Get stand_lon
if (abs(stand_lon) > 180.)
	iter = 0 ;
	while (abs(stand_lon) > 180. && iter < 10)
        if (stand_lon < -180.) 
        	stand_lon = stand_lon + 360.;
        elseif (stand_lon > 180.)
            stand_lon = stand_lon - 360.;
        end
        iter = iter + 1;
	end
end

%Longitude of origin
if ( abs(ref_lon) > 180.)
	iter = 0 ;
	while (abs(ref_lon) > 180. && iter < 10)
        if (ref_lon < -180.)
            ref_lon = ref_lon + 360.;
        elseif (ref_lon > 180.)  
            ref_lon = ref_lon - 360.;
        end
        iter = iter + 1;
	end
end
    
if(truelat1<0)
    hemi=-1.0;
else
    hemi=1.0;
end

if(truelat1~=truelat2)
    k=(log(sind(90-truelat1))-log(sind(90-truelat2)))/...
      (log(tand((90-truelat1)/2))-log(tand((90-truelat2)/2)));
else
    k = sind(abs(truelat1) ) ;
end

delta_lon=ref_lon-stand_lon;
if (delta_lon > 180.)
    delta_lon = delta_lon - 360.;
elseif(delta_lon < -180.)
    delta_lon = delta_lon + 360.;
end

% Calculate mass grid
stagger_i=0;
stagger_j=0;
[lon_m, lat_m,mapfac_m]=WRF_lambert(n,parent_grid_ratio,i_parent_start,j_parent_start,...
                           ref_lat,dx,s_we,e_we,s_sn,e_sn,stand_lon,truelat1,truelat2,...
                           hemi,k,delta_lon,stagger_i,stagger_j);
                       
% Calculate u grid
stagger_i=-0.5;
stagger_j=0;
[lon_u,lat_u,mapfac_u]=WRF_lambert(n,parent_grid_ratio,i_parent_start,j_parent_start,...
                          ref_lat,dx,s_we,e_we,s_sn,e_sn,stand_lon,truelat1,truelat2,...
                          hemi,k,delta_lon,stagger_i,stagger_j);
                   
% Calculate v grid
stagger_i=0;
stagger_j=-0.5;
[lon_v,lat_v,mapfac_v]=WRF_lambert(n,parent_grid_ratio,i_parent_start,j_parent_start,...
                          ref_lat,dx,s_we,e_we,s_sn,e_sn,stand_lon,truelat1,truelat2,...
                          hemi,k,delta_lon,stagger_i,stagger_j);
end
                            
                            
                            
                            
                            
                            
function [lon, lat,mapfac]=WRF_lambert(n,parent_grid_ratio,i_parent_start,j_parent_start,...
                                ref_lat,dx,s_we,e_we,s_sn,e_sn,stand_lon,truelat1,truelat2,...
                                hemi,k,delta_lon,stagger_i,stagger_j)
a=6371000;                % the radius of the earth

M=e_we-s_we+2*abs(stagger_i);              % the nest's full west-east dimension
N=e_sn-s_sn+2*abs(stagger_j);              % the nest's full north-south dimension

rc=a*cosd(truelat1)/k*(tand((90*hemi-ref_lat)/2)/tand((90*hemi-truelat1)/2))^k;
arg=k*delta_lon;
polei=(double(M)+1.)/(2.*double(parent_grid_ratio(n)))-rc*sind(arg)/double(dx);
polej=(double(N)+1.)/(2.*double(parent_grid_ratio(n)))+rc*cosd(arg)/double(dx);
chi1=90. - hemi*truelat1;

M_in_parent = (M - ((real(parent_grid_ratio(n))+1.)/2.))/ real(parent_grid_ratio(n)) + real(i_parent_start(n));
N_in_parent = (N - ((real(parent_grid_ratio(n))+1.)/2.))/ real(parent_grid_ratio(n)) + real(j_parent_start(n));

[i,j]=meshgrid(double(i_parent_start(n)):1/parent_grid_ratio(n):double(M_in_parent),...
               double(j_parent_start(n)):1/parent_grid_ratio(n):double(N_in_parent));
    
I=hemi*i-double(polei);
J=double(polej)-hemi*j;
r=sqrt(double(power(I,2)+power(J,2)));
lon=stand_lon+atand(double(hemi*double(I)./double(J)))/k;
lon=rem(lon+360., 360.);
if(truelat1==truelat2)
    chi=2.0*atand(double(power(( r/(a/dx)/tand(chi1)),(1./k)) * tand((chi1)*0.5)));
else
    chi=2.0*atand(double(power(( r/(a/dx)*k/sind(chi1)),(1./k)) * tand((chi1)*0.5)));
end
lat=(90-chi)*hemi;
lat(r==0)=hemi * 90.;
lon(r==0)=stand_lon;

lon(lon>180)=lon(lon>180)-360;
lon(lon<-180)=lon(lon<-180)+360;

lon=lon';
lat=lat';

% Calculate map factor
theta1         = 90.-truelat1;
theta_m        = 90.-lat;

mapfac       = sind(theta1)./sind(theta_m).*(tand(theta_m/2.)./tand(theta1/2.)).^k;
end