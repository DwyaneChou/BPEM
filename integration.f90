module integration
    use control
    use constant
    use omp_lib
    contains
    subroutine time_integration(ua,va,za,u,v,z,uc,vc,zc,ts,option)
    implicit none
    real    ,intent(in) ::  ua(grid%x_grd_num_u ,grid%y_grd_num_u)  ,&
                            va(grid%x_grd_num_v ,grid%y_grd_num_v)  ,&
                            za(grid%x_grd_num   ,grid%y_grd_num)    ,&
                            u(grid%x_grd_num_u  ,grid%y_grd_num_u)  ,&
                            v(grid%x_grd_num_v  ,grid%y_grd_num_v)  ,&
                            z(grid%x_grd_num    ,grid%y_grd_num)    ,&
                            ts
    integer ,intent(in) ::  option
    real    ,intent(out)::  uc(grid%x_grd_num_u ,grid%y_grd_num_u)     ,&
                            vc(grid%x_grd_num_v ,grid%y_grd_num_v)     ,&
                            zc(grid%x_grd_num   ,grid%y_grd_num)
    
    if(integration_option==1)then
        call time_integration_Square_Conservation(ua,va,za,u,v,z,uc,vc,zc,ts,option)
    elseif(integration_option==2)then
        call time_integration_Energe_Conservation(ua,va,za,u,v,z,uc,vc,zc,ts,option)
    endif
    end subroutine time_integration
    
    !time integrations
    !When option==1, uc=ua+u_t*delta_t
    !When option==2, uc=u_t*delta_t
    !When option==3, uc=u_t
    !ts is the time step
    subroutine time_integration_Square_Conservation(ua,va,za,u,v,z,uc,vc,zc,ts,option)
    implicit none
    real    ,intent(in) ::  ua(grid%x_grd_num_u ,grid%y_grd_num_u)  ,&
                            va(grid%x_grd_num_v ,grid%y_grd_num_v)  ,&
                            za(grid%x_grd_num   ,grid%y_grd_num)    ,&
                            u(grid%x_grd_num_u  ,grid%y_grd_num_u)  ,&
                            v(grid%x_grd_num_v  ,grid%y_grd_num_v)  ,&
                            z(grid%x_grd_num    ,grid%y_grd_num)    ,&
                            ts
    integer ,intent(in) ::  option
    real    ,intent(out)::  uc(grid%x_grd_num_u ,grid%y_grd_num_u)     ,&
                            vc(grid%x_grd_num_v ,grid%y_grd_num_v)     ,&
                            zc(grid%x_grd_num   ,grid%y_grd_num)
    real    zm(grid%x_grd_num     ,grid%y_grd_num)
    real    u_t(grid%x_grd_num_u,grid%y_grd_num_u),&
            v_t(grid%x_grd_num_v,grid%y_grd_num_v),&
            z_t(grid%x_grd_num  ,grid%y_grd_num)
    real    fu(grid%x_grd_num_u,grid%y_grd_num_u)  ,&
            fv(grid%x_grd_num_v,grid%y_grd_num_v)  ,&
            m(grid%x_grd_num    ,grid%y_grd_num)    ,&
            mu(grid%x_grd_num_u,grid%y_grd_num_u)  ,&
            mv(grid%x_grd_num_v,grid%y_grd_num_v)
    real c
    integer i,j
    integer u_xs,&
            u_xe,&
            u_ys,&
            u_ye,&
            v_xs,&
            v_xe,&
            v_ys,&
            v_ye,&
            m_xs,&
            m_xe,&
            m_ys,&
            m_ye

    u_xs    =   2
    u_xe    =   grid%x_grd_num_u-1
    u_ys    =   2
    u_ye    =   grid%y_grd_num_u-1
    v_xs    =   2
    v_xe    =   grid%x_grd_num_v-1
    v_ys    =   2
    v_ye    =   grid%y_grd_num_v-1
    m_xs    =   2
    m_xe    =   grid%x_grd_num-1
    m_ys    =   2
    m_ye    =   grid%y_grd_num-1

    fu  =   grid%f_u
    fv  =   grid%f_v
    m   =   grid%mapfac_m
    mu  =   grid%mapfac_u
    mv  =   grid%mapfac_v
    
    c=0.25d0/grid%dx
    u_t=0.d0
    v_t=0.d0
    z_t=0.d0
    
    zm=z/m
    
    !Loop solve
    !$omp parallel private(i)
    !$omp do
    do j=u_ys,u_ye
        do i=u_xs,u_xe
            u_t(i,j)=-mu(i,j)*c*((u(i+1,j)+u(i,j))*(u(i+1,j)-u(i,j))            &
                                 +(u(i,j)+u(i-1,j))*(u(i,j)-u(i-1,j))           &
                                 +(v(i,j+1)+v(i-1,j+1))*(u(i,j+1)-u(i,j))       &
                                 +(v(i,j)+v(i-1,j))*(u(i,j)-u(i,j-1))           &
                                 +4.d0*g*(z(i,j)-z(i-1,j)))                     &
                                 +0.25d0*(fu(i,j)                               &
                                 +c*((u(i+1,j)+u(i,j))*(mv(i,j+1)-mv(i,j))      &
                                    +(u(i,j)+u(i-1,j))*(mv(i-1,j+1)-mv(i,j))    &
                                    -(v(i,j+1)+v(i,j))*(mu(i+1,j)-mu(i,j))      &
                                    -(v(i-1,j+1)+v(i-1,j))*(mu(i,j)-mu(i-1,j))))&
                                    *(v(i,j+1)+v(i-1,j+1)+v(i,j)+v(i-1,j))
        enddo
    enddo
    !$omp end do
    !$omp do
    do j=v_ys,v_ye
        do i=v_xs,v_xe
            v_t(i,j)=-mv(i,j)*c*((u(i+1,j)+u(i+1,j-1))*(v(i+1,j)-v(i,j))        &
                                 +(u(i,j)+u(i,j-1))*(v(i,j)-v(i-1,j))           &
                                 +(v(i,j+1)+v(i,j))*(v(i,j+1)-v(i,j))           &
                                 +(v(i,j)+v(i,j-1))*(v(i,j)-v(i,j-1))           &
                                 +4.d0*g*(z(i,j)-z(i,j-1)))                     &
                                 -0.25d0*(fv(i,j)                               &
                                 +c*((u(i+1,j)+u(i,j))*(mv(i,j+1)-mv(i,j))      &
                                    +(u(i+1,j-1)+u(i,j-1))*(mv(i,j)-mv(i,j-1))  &
                                    -(v(i,j+1)+v(i,j))*(mu(i+1,j)-mu(i,j))      &
                                    -(v(i,j)+v(i,j-1))*(mu(i+1,j-1)-mu(i,j-1))))&
                                    *(u(i+1,j)+u(i,j)+u(i+1,j-1)+u(i,j-1))
        enddo
    enddo
    !$omp end do
    !$omp do
    do j=m_ys,m_ye
        do i=m_xs,m_xe
            z_t(i,j)=-m(i,j)**2*2.d0*c*(u(i+1,j)*(zm(i+1,j)-zm(i,j))      &
                                       +u(i,j)*(zm(i,j)-zm(i-1,j))        &
                                       +v(i,j+1)*(zm(i,j+1)-zm(i,j))      &
                                       +v(i,j)*(zm(i,j)-zm(i,j-1))        &
                                       +2.d0*zm(i,j)                      &
                                       *(u(i+1,j)-u(i,j)+v(i,j+1)-v(i,j)))
        enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    if(option==1)then
        uc=ua+u_t*ts
        vc=va+v_t*ts
        zc=za+z_t*ts
    elseif(option==2)then
        uc=u_t*ts
        vc=v_t*ts
        zc=z_t*ts
    elseif(option==3)then
        uc=u_t
        vc=v_t
        zc=z_t
    endif

    end subroutine time_integration_Square_Conservation

    
    
    subroutine time_integration_Energe_Conservation(ua,va,za,u,v,z,uc,vc,zc,ts,option)
    implicit none
    real    ,intent(in) ::  ua(grid%x_grd_num_u ,grid%y_grd_num_u)  ,&
                            va(grid%x_grd_num_v ,grid%y_grd_num_v)  ,&
                            za(grid%x_grd_num   ,grid%y_grd_num)    ,&
                            u(grid%x_grd_num_u  ,grid%y_grd_num_u)  ,&
                            v(grid%x_grd_num_v  ,grid%y_grd_num_v)  ,&
                            z(grid%x_grd_num    ,grid%y_grd_num)    ,&
                            ts
    integer ,intent(in) ::  option
    real    ,intent(out)::  uc(grid%x_grd_num_u ,grid%y_grd_num_u)     ,&
                            vc(grid%x_grd_num_v ,grid%y_grd_num_v)     ,&
                            zc(grid%x_grd_num   ,grid%y_grd_num)
    real    zm(grid%x_grd_num     ,grid%y_grd_num)
    real    u_t(grid%x_grd_num_u,grid%y_grd_num_u) ,&
            v_t(grid%x_grd_num_v,grid%y_grd_num_v) ,&
            zu_t(grid%x_grd_num_u,grid%y_grd_num_u),&
            zv_t(grid%x_grd_num_v,grid%y_grd_num_v),&
            z_t(grid%x_grd_num  ,grid%y_grd_num)
    real    fu(grid%x_grd_num_u,grid%y_grd_num_u)  ,&
            fv(grid%x_grd_num_v,grid%y_grd_num_v)  ,&
            m(grid%x_grd_num    ,grid%y_grd_num)    ,&
            mu(grid%x_grd_num_u,grid%y_grd_num_u)  ,&
            mv(grid%x_grd_num_v,grid%y_grd_num_v)
    real c1,c2,c3,c4
    integer i,j
    integer u_xs,&
            u_xe,&
            u_ys,&
            u_ye,&
            v_xs,&
            v_xe,&
            v_ys,&
            v_ye,&
            m_xs,&
            m_xe,&
            m_ys,&
            m_ye

    u_xs    =   3
    u_xe    =   grid%x_grd_num_u-2
    u_ys    =   2
    u_ye    =   grid%y_grd_num_u-1
    v_xs    =   2
    v_xe    =   grid%x_grd_num_v-1
    v_ys    =   3
    v_ye    =   grid%y_grd_num_v-2
    m_xs    =   2
    m_xe    =   grid%x_grd_num-1
    m_ys    =   2
    m_ye    =   grid%y_grd_num-1

    fu  =   grid%f_u
    fv  =   grid%f_v
    m   =   grid%mapfac_m
    mu  =   grid%mapfac_u
    mv  =   grid%mapfac_v
    
    c1=0.125/grid%dx
    c2=g/(2.*grid%dx)
    c3=0.25/grid%dx
    c4=0.5/grid%dx
    
    u_t=0.d0
    v_t=0.d0
    zu_t=0.d0
    zv_t=0.d0
    z_t=0.d0
    
    zm=z/m
    
    !Loop solve
    !$omp parallel private(i)
    !$omp do
    do j=u_ys,u_ye
        do i=u_xs,u_xe
            zu_t(i,j)=-mu(i,j)**2*(c1*(((zm(i+1,j)+zm(i,j))*u(i+1,j)+(zm(i,j)+zm(i-1,j))*u(i,j))*(u(i+1,j)+u(i,j))        &
                                      -((zm(i,j)+zm(i-1,j))*u(i,j)+(zm(i-1,j)+zm(i-2,j))*u(i-1,j))*(u(i,j)+u(i-1,j))      &
                                      +((zm(i,j+1)+zm(i,j))*v(i,j+1)+(zm(i-1,j+1)+zm(i-1,j))*v(i-1,j+1))*(u(i,j+1)+u(i,j))&
                                      -((zm(i,j)+zm(i,j-1))*v(i,j)+(zm(i-1,j)+zm(i-1,j-1))*v(i-1,j))*(u(i,j)+u(i,j-1)))   &
                                    +c2*(zm(i,j)+zm(i-1,j))*(z(i,j)-z(i-1,j)))                                            &
                      +0.125*(fu(i,j)                                                                                     &
                             +c3*((u(i+1,j)+u(i,j))*(mv(i,j+1)-mv(i,j))                                                   &
                                 +(u(i,j)+u(i-1,j))*(mv(i-1,j+1)-mv(i-1,j))                                               &
                                 -(v(i,j+1)+v(i,j))*(mu(i+1,j)-mu(i,j))                                                   &
                                 -(v(i-1,j+1)+v(i-1,j))*(mu(i,j)-mu(i-1,j))))                                             &
                                *((z(i,j+1)+z(i,j))*v(i,j+1)                                                              &
                                 +(z(i-1,j+1)+z(i-1,j))*v(i-1,j+1)                                                        &
                                 +(z(i,j)+z(i,j-1))*v(i,j)                                                                &
                                 +(z(i-1,j)+z(i-1,j-1))*v(i-1,j))
        enddo
    enddo
    !$omp end do
    !$omp do
    do j=v_ys,v_ye
        do i=v_xs,v_xe
            zv_t(i,j)=-mv(i,j)**2*(c1*(((zm(i+1,j)+zm(i,j))*u(i+1,j)+(zm(i+1,j-1)+zm(i,j-1))*u(i+1,j-1))*(v(i+1,j)+v(i,j))&
                                      -((zm(i,j)+zm(i-1,j))*u(i,j)+(zm(i,j-1)+zm(i-1,j-1))*u(i,j-1))*(v(i,j)+v(i-1,j))    &
                                      +((zm(i,j+1)+zm(i,j))*v(i,j+1)+(zm(i,j)+zm(i,j-1))*v(i,j))*(v(i,j+1)+v(i,j))        &
                                      -((zm(i,j)+zm(i,j-1))*v(i,j)+(zm(i,j-1)+zm(i,j-2))*v(i,j-1))*(v(i,j)+v(i,j-1)))     &
                                    +c2*(zm(i,j)+zm(i,j-1))*(z(i,j)-z(i,j-1)))                                            &
                      -0.125*(fv(i,j)                                                                                     &
                             +c3*((u(i+1,j)+u(i,j))*(mv(i,j+1)-mv(i,j))                                                   &
                                 +(u(i+1,j-1)+u(i,j-1))*(mv(i,j)-mv(i,j-1))                                               &
                                 -(v(i,j+1)+v(i,j))*(mu(i+1,j)-mu(i,j))                                                   &
                                 -(v(i,j)+v(i,j-1))*(mu(i+1,j-1)-mu(i,j-1))))                                             &
                                *((z(i+1,j)+z(i,j))*u(i+1,j)                                                              &
                                 +(z(i+1,j-1)+z(i,j-1))*u(i+1,j-1)                                                        &
                                 +(z(i,j)+z(i-1,j))*u(i,j)                                                                &
                                 +(z(i,j-1)+z(i-1,j-1))*u(i,j-1))
        enddo
    enddo
    !$omp end do
    !$omp do
    do j=m_ys,m_ye
        do i=m_xs,m_xe
            z_t(i,j)=-m(i,j)**2*c4*((zm(i+1,j)+zm(i,j))*u(i+1,j)&
                                   -(zm(i,j)+zm(i-1,j))*u(i,j)  &
                                   +(zm(i,j+1)+zm(i,j))*v(i,j+1)&
                                   -(zm(i,j)+zm(i,j-1))*v(i,j))
        enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    
    !$omp parallel private(i)
    !$omp do
    do j=u_ys,u_ye
        do i=u_xs,u_xe
            u_t(i,j)=2./(z(i,j)+z(i-1,j))*(zu_t(i,j)-u(i,j)*(z_t(i,j)+z_t(i-1,j))*0.5)
        enddo
    enddo
    !$omp end do
    !$omp do
    do j=v_ys,v_ye
        do i=v_xs,v_xe
            v_t(i,j)=2./(z(i,j)+z(i,j-1))*(zv_t(i,j)-v(i,j)*(z_t(i,j)+z_t(i,j-1))*0.5)
        enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    !print*,'max u_t = ',maxval(abs(u_t)),' max v_t = ',maxval(abs(v_t)),' max z_t = ',maxval(abs(z_t)),' maxloc=',maxloc(z_t)
    !pause
    
    if(option==1)then
        uc=ua+u_t*ts
        vc=va+v_t*ts
        zc=za+z_t*ts
    elseif(option==2)then
        uc=u_t*ts
        vc=v_t*ts
        zc=z_t*ts
    elseif(option==3)then
        uc=u_t
        vc=v_t
        zc=z_t
    endif

    end subroutine time_integration_Energe_Conservation
    
end module integration