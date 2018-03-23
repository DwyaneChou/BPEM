clc
clear
BPEM_output='E:\study\BPEM\output\BPEM_output.nc';
picture_output_path='E:\study\BPEM\postprocess\pictures';
obs_file='E:\study\BPEM\initialize\BPEM_input.nc';

plot_obs=1;%1 for true 0 for false

map_line_width=0.1;
adjust_subplot_x_space=0.037;
title_font_size=6;
adjust_title_y_position=0.15;
adjust_color_bar_y_position=1;

system='Windows'; % Choose from Windows or Linux
if strcmp(system,'Windows')
    slash='\';
elseif strcmp(system,'Linux')
    slash='/';
end

XLONG_M = ncread(BPEM_output,'XLONG_M');
XLAT_M  = ncread(BPEM_output,'XLAT_M');
XLONG_U = ncread(BPEM_output,'XLONG_U');
XLAT_U  = ncread(BPEM_output,'XLAT_U');
XLONG_V = ncread(BPEM_output,'XLONG_V');
XLAT_V  = ncread(BPEM_output,'XLAT_V');
Times   = ncread(BPEM_output,'Times');
Z       = ncread(BPEM_output,'Z');
U       = ncread(BPEM_output,'U');
V       = ncread(BPEM_output,'V');

CEN_LON = ncreadatt(BPEM_output,'/','CEN_LON');
CEN_LAT = ncreadatt(BPEM_output,'/','CEN_LAT');

XLONG_M(XLONG_M<0)=360+XLONG_M(XLONG_M<0);
XLONG_U(XLONG_U<0)=360+XLONG_U(XLONG_U<0);
XLONG_V(XLONG_V<0)=360+XLONG_V(XLONG_V<0);

time_length=size(Times,2);

if plot_obs
    U_obs=ncread(obs_file,'U');
    V_obs=ncread(obs_file,'V');
    Z_obs=ncread(obs_file,'Z');
end

%%%%%%%%
% Plot %
%%%%%%%%
ChinaL=shaperead('D:\MeteoInfo\map\bou2_4p.shp');
bou2_4px=[ChinaL(:).X];
bou2_4py=[ChinaL(:).Y];
parfor i=1:time_length
    figure('visible','off')
    subplot(1,2,1)
    ax          = gca;
    axpos       = ax.Position;
    axpos(1)    = axpos(1)+adjust_subplot_x_space;
    ax.Position = axpos;
    
    plot_X_interval   = max(CEN_LON-min(min(XLONG_U(:,:,i))),max(max(XLONG_U(:,:,i)))-CEN_LON);
    plot_Y_interval   = max(CEN_LAT-min(min(XLAT_U(:,:,i))),max(max(XLAT_U(:,:,i)))-CEN_LAT);
    plot_right_lon    = double(CEN_LON+plot_X_interval);
    plot_left_lon     = double(CEN_LON-plot_X_interval);
    plot_upper_lat    = double(CEN_LAT+plot_Y_interval);
    plot_bottom_lat   = double(CEN_LAT-plot_Y_interval);
    m_proj('lambert','lon',[plot_left_lon,plot_right_lon],'lat',[plot_bottom_lat,plot_upper_lat]);
    m_pcolor(double(XLONG_U(:,:,i)),double(XLAT_U(:,:,i)),double(U(:,:,i)));
    hold on
    m_plot(bou2_4px,bou2_4py,'k','linewidth',map_line_width);
    m_grid('fontsize',4);
    m_coast('color','k','linewidth',map_line_width);
    colormap(jet)
    shading interp
    set(gca,'fontsize',5,'CLim',[-30,50])
    
%   Plot colorbar and change colorbar width
    c           = colorbar('southoutside','Ticks',-100:5:100);
    cpos        = c.Position;
    cpos(2)     = cpos(2)-adjust_color_bar_y_position*cpos(4);
    cpos(4)     = 0.5*cpos(4);
    c.Position  = cpos;
    
%   Plot title and change title width
    t=title(['U at 500hPa ',Times(1:10,i)',' ',Times(12:end,i)'],'fontsize',title_font_size,'fontweight','bold');
    tpos        = t.Position;
    tpos(2)     = tpos(2)+adjust_title_y_position;
    t.Position  = tpos;
    
    if plot_obs
        subplot(1,2,2)
        ax          = gca;
        axpos       = ax.Position;
        axpos(1)    = axpos(1)-adjust_subplot_x_space;
        ax.Position = axpos;
        
        plot_X_interval   = max(CEN_LON-min(min(XLONG_U(:,:,i))),max(max(XLONG_U(:,:,i)))-CEN_LON);
        plot_Y_interval   = max(CEN_LAT-min(min(XLAT_U(:,:,i))),max(max(XLAT_U(:,:,i)))-CEN_LAT);
        plot_right_lon    = double(CEN_LON+plot_X_interval);
        plot_left_lon     = double(CEN_LON-plot_X_interval);
        plot_upper_lat    = double(CEN_LAT+plot_Y_interval);
        plot_bottom_lat   = double(CEN_LAT-plot_Y_interval);
        m_proj('lambert','lon',[plot_left_lon,plot_right_lon],'lat',[plot_bottom_lat,plot_upper_lat]);
        m_pcolor(double(XLONG_U(:,:,i)),double(XLAT_U(:,:,i)),double(U_obs(:,:,i)));
        hold on
        m_plot(bou2_4px,bou2_4py,'k','linewidth',map_line_width);
        m_grid('fontsize',4);
        m_coast('color','k','linewidth',map_line_width);
        colormap(jet)
        shading interp
        set(gca,'fontsize',5,'CLim',[-30,50])
        
        %   Plot colorbar and change colorbar width
        c           = colorbar('southoutside','Ticks',-100:5:100);
        cpos        = c.Position;
        cpos(2)     = cpos(2)-adjust_color_bar_y_position*cpos(4);
        cpos(4)     = 0.5*cpos(4);
        c.Position  = cpos;
        
        %   Plot title and change title width
        t=title(['U obs at 500hPa ',Times(1:10,i)',' ',Times(12:end,i)'],'fontsize',title_font_size,'fontweight','bold');
        tpos        = t.Position;
        tpos(2)     = tpos(2)+adjust_title_y_position;
        t.Position  = tpos;
    end
    
    print('-dpng','-zbuffer','-r300',[picture_output_path,slash,'u',num2str(i,'%04d'),'.png']);
end