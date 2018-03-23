module subroutines
    use control
    use constant
    contains
    !space smoothing for internal points 区域内5点平滑(正逆平滑)
    !l=1为只执行正平滑，l=2为执行正逆平滑
    subroutine smooth_5_points(a,smooth_coefficient,x_grd_num,y_grd_num,l)
    implicit none
    integer ,intent(in)     ::  x_grd_num,y_grd_num
    real    ,intent(in)     ::  smooth_coefficient
    real    ,intent(inout)  ::  a(x_grd_num,y_grd_num)
    real                        w(x_grd_num,y_grd_num)
    integer                     l,xs,xe,ys,ye
    integer                     i,j

    xs=2
    ys=2
    xe=x_grd_num-1
    ye=y_grd_num-1

    do i=xs,xe
        do j=xs,ye
            w(i,j)=a(i,j)+smooth_coefficient*0.25*(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1)-4.*a(i,j))
        enddo
    enddo
    a(xs:xe,ys:ye)=w(xs:xe,ys:ye)

    if(l==2) then
        do i=xs,xe
            do j=ys,ye
                w(i,j)=a(i,j)-smooth_coefficient*0.25*(a(i-1,j)+a(i+1,j)+a(i,j-1)+a(i,j+1)-4.*a(i,j))
            enddo
        enddo
        a(xs:xe,ys:ye)=w(xs:xe,ys:ye)
    end if
    end subroutine smooth_5_points

    !transmiting arrays  数组传送
    subroutine transmit_arrays(ua,va,za,ub,vb,zb)
    real    ,intent(in) ::  ub(grid%x_grd_num_u ,grid%y_grd_num_u),&
                            vb(grid%x_grd_num_v ,grid%y_grd_num_v),&
                            zb(grid%x_grd_num   ,grid%y_grd_num)
    real    ,intent(out)::  ua(grid%x_grd_num_u  ,grid%y_grd_num_u),&
                            va(grid%x_grd_num_v  ,grid%y_grd_num_v),&
                            za(grid%x_grd_num    ,grid%y_grd_num)
    ua=ub
    va=vb
    za=zb
    end subroutine transmit_arrays

    !Lateral Boundary Conditions by Davies and Turner (1977) modified by Skamarock(2008)
    subroutine WRF_LBC(var_in,var_t,var_BYS,var_BYE,var_BXS,var_BXE,var_BTYS,var_BTYE,var_BTXS,var_BTXE ,&
                      F1,F2,interval_seconds,integrated_step_num,ts,x_grd_num,y_grd_num                 ,&
                      bdy_width,spec_zone,bdy_time_num,bdy_tend_time_num)
    implicit none
    integer,intent(in)::spec_zone           ,&
                        x_grd_num           ,&
                        y_grd_num           ,&
                        bdy_width           ,&
                        bdy_time_num        ,&
                        bdy_tend_time_num
    
    real   ,intent(in)::var_in (x_grd_num,y_grd_num)                   ,&
                        var_BYS(x_grd_num,bdy_width,bdy_time_num)      ,&
                        var_BYE(x_grd_num,bdy_width,bdy_time_num)      ,&
                        var_BXS(y_grd_num,bdy_width,bdy_time_num)      ,&
                        var_BXE(y_grd_num,bdy_width,bdy_time_num)      ,&
                        var_BTYS(x_grd_num,bdy_width,bdy_tend_time_num),&
                        var_BTYE(x_grd_num,bdy_width,bdy_tend_time_num),&
                        var_BTXS(y_grd_num,bdy_width,bdy_tend_time_num),&
                        var_BTXE(y_grd_num,bdy_width,bdy_tend_time_num),&
                        F1(spec_zone+1:bdy_width-1)                    ,&
                        F2(spec_zone+1:bdy_width-1)                    ,&
                        interval_seconds                               ,&
                        integrated_step_num                            ,&
                        ts

    real  ,intent(out)::var_t(x_grd_num,y_grd_num)

    real                BYS        (x_grd_num,bdy_width)    ,&
                        BYE        (x_grd_num,bdy_width)    ,&
                        BXS        (y_grd_num,bdy_width)    ,&
                        BXE        (y_grd_num,bdy_width)    ,&
                        dBYS       (x_grd_num,bdy_width)    ,&
                        dBYE       (x_grd_num,bdy_width)    ,&
                        dBXS       (y_grd_num,bdy_width)    ,&
                        dBXE       (y_grd_num,bdy_width)

    integer             i,j,k,i_inv,j_inv,n
    integer             bdy_stage,bdy_stage_index
    integer             jm1,jp1
    real                ij_angle
    
    var_t       = 0.
    
    bdy_stage       = int(integrated_step_num*time_step/interval_seconds)+1
    bdy_stage_index = min(grid%BDY_TEND_TIME_NUM,bdy_stage)
    BYS = var_BYS(:,:,bdy_stage)+dble(integrated_step_num*time_step-(bdy_stage-1)*interval_seconds)*var_BTYS(:,:,bdy_stage_index)/3600.
    BYE = var_BYE(:,:,bdy_stage)+dble(integrated_step_num*time_step-(bdy_stage-1)*interval_seconds)*var_BTYE(:,:,bdy_stage_index)/3600.
    BXS = var_BXS(:,:,bdy_stage)+dble(integrated_step_num*time_step-(bdy_stage-1)*interval_seconds)*var_BTXS(:,:,bdy_stage_index)/3600.
    BXE = var_BXE(:,:,bdy_stage)+dble(integrated_step_num*time_step-(bdy_stage-1)*interval_seconds)*var_BTXE(:,:,bdy_stage_index)/3600.
    
    forall(i=1:x_grd_num,j=1:bdy_width)
        dBYS(i,j) = BYS(i,j)-var_in(i,j)
        dBYE(i,j) = BYE(i,j)-var_in(i,y_grd_num-j+1)
    endforall
    
    forall(i=1:bdy_width,j=1:y_grd_num)
        dBXS(j,i) = BXS(j,i)-var_in(i,j)
        dBXE(j,i) = BXE(j,i)-var_in(x_grd_num-i+1,j)
    endforall
    
    do j=1,x_grd_num
        jm1 = max(j-1,1)
        jp1 = min(j+1,x_grd_num)
        forall(i=spec_zone+1:bdy_width-spec_zone)
            var_t(j,i)              =   F1(i)/ts*dBYS(j,i)-F2(i)/ts*(dBYS(jm1,i)+dBYS(jp1,i)+dBYS(j,i+1)+dBYS(j,i-1)-4.*dBYS(j,i))
            var_t(j,y_grd_num-i+1)  =   F1(i)/ts*dBYE(j,i)-F2(i)/ts*(dBYE(jm1,i)+dBYE(jp1,i)+dBYE(j,i+1)+dBYE(j,i-1)-4.*dBYE(j,i))
        endforall
    enddo
    do j=1,y_grd_num
        jm1 = max(j-1,1)
        jp1 = min(j+1,y_grd_num)
        forall(i=spec_zone+1:bdy_width-spec_zone)
            var_t(i,j)              =   F1(i)/ts*dBXS(j,i)-F2(i)/ts*(dBXS(jm1,i)+dBXS(jp1,i)+dBXS(j,i+1)+dBXS(j,i-1)-4.*dBXS(j,i))
            var_t(x_grd_num-i+1,j)  =   F1(i)/ts*dBXE(j,i)-F2(i)/ts*(dBXE(jm1,i)+dBXE(jp1,i)+dBXE(j,i+1)+dBXE(j,i-1)-4.*dBXE(j,i))
        endforall
    enddo
    
    do i=spec_zone+1,bdy_width-1
        do j=spec_zone+1,bdy_width-1
            ij_angle=atan(dble(j)/dble(i))
            i_inv=x_grd_num-i+1
            j_inv=y_grd_num-j+1
            var_t(i,j)          =   (sin(ij_angle)**2*F1(i)+cos(ij_angle)**2*F1(j))/ts*dBYS(i,j)    -(sin(ij_angle)**2*F2(i)+cos(ij_angle)**2*F2(j))/ts*(dBYS(i-1,j)    +dBYS(i+1,j)    +dBYS(i,j+1)    +dBYS(i,j-1)    -4.*dBYS(i,j))
            var_t(i_inv,j)      =   (sin(ij_angle)**2*F1(i)+cos(ij_angle)**2*F1(j))/ts*dBXE(j,i)    -(sin(ij_angle)**2*F2(i)+cos(ij_angle)**2*F2(j))/ts*(dBXE(j-1,i)    +dBXE(j+1,i)    +dBXE(j,i+1)    +dBXE(j,i-1)    -4.*dBXE(j,i))
            var_t(i,j_inv)      =   (sin(ij_angle)**2*F1(i)+cos(ij_angle)**2*F1(j))/ts*dBYE(i,j)    -(sin(ij_angle)**2*F2(i)+cos(ij_angle)**2*F2(j))/ts*(dBYE(i-1,j)    +dBYE(i+1,j)    +dBYE(i,j+1)    +dBYE(i,j-1)    -4.*dBYE(i,j))
            var_t(i_inv,j_inv)  =   (sin(ij_angle)**2*F1(i)+cos(ij_angle)**2*F1(j))/ts*dBYE(i_inv,j)-(sin(ij_angle)**2*F2(i)+cos(ij_angle)**2*F2(j))/ts*(dBYE(i_inv-1,j)+dBYE(i_inv+1,j)+dBYE(i_inv,j+1)+dBYE(i_inv,j-1)-4.*dBYE(i_inv,j))
        enddo
    enddo
    
    forall (i=1:spec_zone,j=1:x_grd_num)
        var_t(j,i)              =   dBYS(j,i)/ts
        var_t(j,y_grd_num-i+1)  =   dBYE(j,i)/ts
    endforall
    forall (i=1:spec_zone,j=1:y_grd_num)
        var_t(i,j)              =   dBXS(j,i)/ts
        var_t(x_grd_num-i+1,j)  =   dBXE(j,i)/ts
    endforall
    
    end subroutine WRF_LBC
    
    ! Calculate Lateral Boundary Condition Weight Function Coefficients
    subroutine LBC_WFC(F1,F2,bdy_width,spec_zone,RelaxZone,dx)
        implicit none
        integer,intent(in)  ::  bdy_width,spec_zone,RelaxZone
        real,intent(in)     ::  dx
        real,intent(out)    ::  F1(spec_zone+1:bdy_width-1),&
                                F2(spec_zone+1:bdy_width-1)
        integer n
        
        forall (n=spec_zone+1:bdy_width-1)
            F1(n)= 0.1d0*dble(bdy_width-n)/dble(RelaxZone-1.d0)*exp(-dble(n-spec_zone-1.d0)*spec_exp) !spec_exp is a parameter defined by user, 0 for linear ramp, 0.33 =~ 3*dx exp decay factor
            F2(n)=0.02d0*dble(bdy_width-n)/dble(RelaxZone-1.d0)*exp(-dble(n-spec_zone-1.d0)*spec_exp) !spec_exp is a parameter defined by user, 0 for linear ramp, 0.33 =~ 3*dx exp decay factor
        endforall
    end subroutine LBC_WFC

    
    
    
    !Update Lateral Boundary Conditions by using WRF scheme
    subroutine ULBC_WRF(u,v,z,integerated_step_num,ts)
    implicit none
    real    ,intent(inout)::u(grid%x_grd_num_u,grid%y_grd_num_u),&
                            v(grid%x_grd_num_v,grid%y_grd_num_v),&
                            z(grid%x_grd_num  ,grid%y_grd_num  )
    
    real    ,intent(in)::   integerated_step_num
    real    ,intent(in)::   ts !time_step used by updating lbc
    
    real                    u_bdy_temp(grid%x_grd_num_u,grid%y_grd_num_u),&
                            v_bdy_temp(grid%x_grd_num_v,grid%y_grd_num_v),&
                            z_bdy_temp(grid%x_grd_num  ,grid%y_grd_num  )
    integer                 iubs,iube,jubs,jube,iuts,iute,juts,jute,iuls,iule,juls,jule,iurs,iure,jurs,jure,&
                            ivbs,ivbe,jvbs,jvbe,ivts,ivte,jvts,jvte,ivls,ivle,jvls,jvle,ivrs,ivre,jvrs,jvre,&
                            izbs,izbe,jzbs,jzbe,izts,izte,jzts,jzte,izls,izle,jzls,jzle,izrs,izre,jzrs,jzre
    
    !$omp parallel
    !$omp sections    
    !$omp section
    !u bottom bdy
    iubs=1
    iube=grid%x_grd_num_u
    jubs=1
    jube=grid%bdy_width
    !$omp section
    !u top bdy
    iuts=1
    iute=grid%x_grd_num_u
    juts=grid%y_grd_num_u-grid%bdy_width+1
    jute=grid%y_grd_num_u
    !$omp section
    !u left bdy
    iuls=1
    iule=grid%bdy_width
    juls=1
    jule=grid%y_grd_num_u
    !$omp section
    !u right bdy
    iurs=grid%x_grd_num_u-grid%bdy_width+1
    iure=grid%x_grd_num_u
    jurs=1
    jure=grid%y_grd_num_u
    !$omp section
    !v bottom bdy
    ivbs=1
    ivbe=grid%x_grd_num_v
    jvbs=1
    jvbe=grid%bdy_width
    !$omp section
    !v top bdy
    ivts=1
    ivte=grid%x_grd_num_v
    jvts=grid%y_grd_num_v-grid%bdy_width+1
    jvte=grid%y_grd_num_v
    !$omp section
    !v left bdy
    ivls=1
    ivle=grid%bdy_width
    jvls=1
    jvle=grid%y_grd_num_v
    !$omp section
    !v right bdy
    ivrs=grid%x_grd_num_v-grid%bdy_width+1
    ivre=grid%x_grd_num_v
    jvrs=1
    jvre=grid%y_grd_num_v
    !$omp section
    !z bottom bdy
    izbs=1
    izbe=grid%x_grd_num
    jzbs=1
    jzbe=grid%bdy_width
    !$omp section
    !z top bdy
    izts=1
    izte=grid%x_grd_num
    jzts=grid%y_grd_num-grid%bdy_width+1
    jzte=grid%y_grd_num
    !$omp section
    !z left bdy
    izls=1
    izle=grid%bdy_width
    jzls=1
    jzle=grid%y_grd_num
    !$omp section
    !z right bdy
    izrs=grid%x_grd_num-grid%bdy_width+1
    izre=grid%x_grd_num
    jzrs=1
    jzre=grid%y_grd_num
    !$omp end sections
    !$omp end parallel
    
    !$omp parallel
    !$omp sections
    !$omp section
    call WRF_LBC(u,u_bdy_temp,grid%U_BYS,grid%U_BYE,grid%U_BXS,grid%U_BXE,grid%U_BTYS,grid%U_BTYE,grid%U_BTXS,grid%U_BTXE  ,&
                grid%F1,grid%F2,grid%interval_seconds,integerated_step_num,ts,grid%x_grd_num_u,grid%y_grd_num_u            ,&
                grid%bdy_width,grid%Spec_Zone,grid%bdy_time_num,grid%bdy_tend_time_num)
    !$omp section
    call WRF_LBC(v,v_bdy_temp,grid%V_BYS,grid%V_BYE,grid%V_BXS,grid%V_BXE,grid%V_BTYS,grid%V_BTYE,grid%V_BTXS,grid%V_BTXE  ,&
                grid%F1,grid%F2,grid%interval_seconds,integerated_step_num,ts,grid%x_grd_num_v,grid%y_grd_num_v            ,&
                grid%bdy_width,grid%Spec_Zone,grid%bdy_time_num,grid%bdy_tend_time_num)
    !$omp section
    call WRF_LBC(z,z_bdy_temp,grid%Z_BYS,grid%Z_BYE,grid%Z_BXS,grid%Z_BXE,grid%Z_BTYS,grid%Z_BTYE,grid%Z_BTXS,grid%Z_BTXE  ,&
                grid%F1,grid%F2,grid%interval_seconds,integerated_step_num,ts,grid%x_grd_num,grid%y_grd_num                ,&
                grid%bdy_width,grid%spec_zone,grid%bdy_time_num,grid%bdy_tend_time_num)
    !$omp end sections
    !$omp end parallel
    
    !$omp parallel
    !$omp sections    
    !$omp section
    !u bottom bound
    u(iubs:iube,jubs:jube)=u(iubs:iube,jubs:jube)+ts*u_bdy_temp(iubs:iube,jubs:jube)
    !$omp section
    !u top bound
    u(iuts:iute,juts:jute)=u(iuts:iute,juts:jute)+ts*u_bdy_temp(iuts:iute,juts:jute)
    !$omp section
    !u left bound
    u(iuls:iule,juls:jule)=u(iuls:iule,juls:jule)+ts*u_bdy_temp(iuls:iule,juls:jule)
    !$omp section
    !u right bound
    u(iurs:iure,jurs:jure)=u(iurs:iure,jurs:jure)+ts*u_bdy_temp(iurs:iure,jurs:jure)
    !$omp section
    !v bottom bound
    v(ivbs:ivbe,jvbs:jvbe)=v(ivbs:ivbe,jvbs:jvbe)+ts*v_bdy_temp(ivbs:ivbe,jvbs:jvbe)
    !$omp section
    !v top bound
    v(ivts:ivte,jvts:jvte)=v(ivts:ivte,jvts:jvte)+ts*v_bdy_temp(ivts:ivte,jvts:jvte)
    !$omp section
    !v left bound
    v(ivls:ivle,jvls:jvle)=v(ivls:ivle,jvls:jvle)+ts*v_bdy_temp(ivls:ivle,jvls:jvle)
    !$omp section
    !v right bound
    v(ivrs:ivre,jvrs:jvre)=v(ivrs:ivre,jvrs:jvre)+ts*v_bdy_temp(ivrs:ivre,jvrs:jvre)
    !$omp section
    !z bottom bound
    z(izbs:izbe,jzbs:jzbe)=z(izbs:izbe,jzbs:jzbe)+ts*z_bdy_temp(izbs:izbe,jzbs:jzbe)
    !$omp section
    !z top bound
    z(izts:izte,jzts:jzte)=z(izts:izte,jzts:jzte)+ts*z_bdy_temp(izts:izte,jzts:jzte)
    !$omp section
    !z left bound
    z(izls:izle,jzls:jzle)=z(izls:izle,jzls:jzle)+ts*z_bdy_temp(izls:izle,jzls:jzle)
    !$omp section
    !z right bound
    z(izrs:izre,jzrs:jzre)=z(izrs:izre,jzrs:jzre)+ts*z_bdy_temp(izrs:izre,jzrs:jzre)
    !$omp end sections
    !$omp end parallel
    
    end subroutine ULBC_WRF
    
    subroutine diffusion(var,x_grd_num,y_grd_num)
        implicit none
        !ts is the time step
        integer ,intent(in)::       x_grd_num,y_grd_num
        real    ,intent(inout)::    var(x_grd_num,y_grd_num)
        real                        diffusion_term(x_grd_num,y_grd_num)
        real                        k2
        integer                     i,j
        
        diffusion_term=0.d0
        !2nd order
        forall(i=2:2,j=2:y_grd_num-1)
            diffusion_term(i,j)=0.125*(var(i+1,j)+var(i-1,j)+var(i,j+1)+var(i,j-1)-4*var(i,j))
        endforall
        
        forall(i=grid%bdy_width+1:x_grd_num-1,j=2:2)
            diffusion_term(i,j)=0.125*(var(i+1,j)+var(i-1,j)+var(i,j+1)+var(i,j-1)-4*var(i,j))
        endforall
        
        forall(i=x_grd_num-1:x_grd_num-1,j=grid%bdy_width+1:y_grd_num-1)
            diffusion_term(i,j)=0.125*(var(i+1,j)+var(i-1,j)+var(i,j+1)+var(i,j-1)-4*var(i,j))
        endforall
        
        forall(i=grid%bdy_width+1:x_grd_num-grid%bdy_width,j=y_grd_num-1:y_grd_num-1)
            diffusion_term(i,j)=0.125*(var(i+1,j)+var(i-1,j)+var(i,j+1)+var(i,j-1)-4*var(i,j))
        endforall
        
        !4th order
        forall(i=3:3,j=3:y_grd_num-2)
            diffusion_term(i,j)=-0.03125*(var(i+2,j)+var(i-2,j)+var(i,j+2)+var(i,j-2)-4*(var(i+1,j)+var(i-1,j)+var(i,j+1)+var(i,j-1))+12*var(i,j))
        endforall
        
        forall(i=3:x_grd_num-2,j=3:3)
            diffusion_term(i,j)=-0.03125*(var(i+2,j)+var(i-2,j)+var(i,j+2)+var(i,j-2)-4*(var(i+1,j)+var(i-1,j)+var(i,j+1)+var(i,j-1))+12*var(i,j))
        endforall
        
        forall(i=x_grd_num-2:x_grd_num-2,j=3:y_grd_num-2)
            diffusion_term(i,j)=-0.03125*(var(i+2,j)+var(i-2,j)+var(i,j+2)+var(i,j-2)-4*(var(i+1,j)+var(i-1,j)+var(i,j+1)+var(i,j-1))+12*var(i,j))
        endforall
        
        forall(i=3:x_grd_num-2,j=y_grd_num-2:y_grd_num-2)
            diffusion_term(i,j)=-0.03125*(var(i+2,j)+var(i-2,j)+var(i,j+2)+var(i,j-2)-4*(var(i+1,j)+var(i-1,j)+var(i,j+1)+var(i,j-1))+12*var(i,j))
        endforall
        
        !6th order
        forall(i=4:grid%bdy_width,j=4:y_grd_num-3)
            diffusion_term(i,j)=0.0078125*(var(i+3,j)+var(i-3,j)+var(i,j+3)+var(i,j-3)-6*(var(i+2,j)+var(i-2,j)+var(i,j+2)+var(i,j-2))+15*(var(i+1,j)+var(i-1,j)+var(i,j+1)+var(i,j-1))-40*var(i,j))
        endforall
        
        forall(i=grid%bdy_width+1:x_grd_num-3,j=4:grid%bdy_width)
            diffusion_term(i,j)=0.0078125*(var(i+3,j)+var(i-3,j)+var(i,j+3)+var(i,j-3)-6*(var(i+2,j)+var(i-2,j)+var(i,j+2)+var(i,j-2))+15*(var(i+1,j)+var(i-1,j)+var(i,j+1)+var(i,j-1))-40*var(i,j))
        endforall
        
        forall(i=x_grd_num-grid%bdy_width+1:x_grd_num-3,j=grid%bdy_width+1:y_grd_num-3)
            diffusion_term(i,j)=0.0078125*(var(i+3,j)+var(i-3,j)+var(i,j+3)+var(i,j-3)-6*(var(i+2,j)+var(i-2,j)+var(i,j+2)+var(i,j-2))+15*(var(i+1,j)+var(i-1,j)+var(i,j+1)+var(i,j-1))-40*var(i,j))
        endforall
        
        forall(i=grid%bdy_width+1:x_grd_num-grid%bdy_width,j=y_grd_num-grid%bdy_width+1:y_grd_num-3)
            diffusion_term(i,j)=0.0078125*(var(i+3,j)+var(i-3,j)+var(i,j+3)+var(i,j-3)-6*(var(i+2,j)+var(i-2,j)+var(i,j+2)+var(i,j-2))+15*(var(i+1,j)+var(i-1,j)+var(i,j+1)+var(i,j-1))-40*var(i,j))
        endforall
        
        !6th order diffusion full field
        !forall(i=4:x_grd_num-3,j=4:y_grd_num-3)
        !    diffusion_term(i,j)=-0.0078125*(var(i+3,j)+var(i-3,j)+var(i,j+3)+var(i,j-3)-6*(var(i+2,j)+var(i-2,j)+var(i,j+2)+var(i,j-2))+15*(var(i+1,j)+var(i-1,j)+var(i,j+1)+var(i,j-1))-40*var(i,j))
        !endforall
        
        !8th order diffusion full field
        forall(i=5:x_grd_num-4,j=5:y_grd_num-4)
            diffusion_term(i,j)=-0.001953125*(var(i+4,j)+var(i-4,j)+var(i,j+4)+var(i,j-4)-8*(var(i+3,j)+var(i-3,j)+var(i,j+3)+var(i,j-3))+28*(var(i+2,j)+var(i-2,j)+var(i,j+2)+var(i,j-2))-56*(var(i+1,j)+var(i-1,j)+var(i,j+1)+var(i,j-1))+140*var(i,j))
        endforall
        
        var=var+diffusion_term
        
    end subroutine diffusion
    
end module subroutines