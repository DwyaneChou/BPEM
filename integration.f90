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
        call time_integration_loop(ua,va,za,u,v,z,uc,vc,zc,ts,option)
    elseif(integration_option==2)then
        call time_integration_matrix(ua,va,za,u,v,z,uc,vc,zc,ts,option)
    endif
    end subroutine time_integration
    
    !time integrations
    !When option==1, uc=ua+u_t*delta_t
    !When option==2, uc=u_t*delta_t
    !When option==3, uc=u_t
    !ts is the time step
    subroutine time_integration_matrix(ua,va,za,u,v,z,uc,vc,zc,ts,option)
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
    
    !Define Matrix Variables
    real    uip1    (grid%x_grd_num_u   ,grid%y_grd_num_u)  ,&
            uim1    (grid%x_grd_num_u   ,grid%y_grd_num_u)  ,&
            ujp1    (grid%x_grd_num_u   ,grid%y_grd_num_u)  ,&
            ujm1    (grid%x_grd_num_u   ,grid%y_grd_num_u)  ,&
            uip1jp1 (grid%x_grd_num_u   ,grid%y_grd_num_u)  ,&
            uip1jm1 (grid%x_grd_num_u   ,grid%y_grd_num_u)  ,&
            uim1jp1 (grid%x_grd_num_u   ,grid%y_grd_num_u)  ,&
            uim1jm1 (grid%x_grd_num_u   ,grid%y_grd_num_u)  ,&
            muip1   (grid%x_grd_num_u   ,grid%y_grd_num_u)  ,&
            muim1   (grid%x_grd_num_u   ,grid%y_grd_num_u)  ,&
            mujp1   (grid%x_grd_num_u   ,grid%y_grd_num_u)  ,&
            mujm1   (grid%x_grd_num_u   ,grid%y_grd_num_u)  ,&
            muip1jp1(grid%x_grd_num_u   ,grid%y_grd_num_u)  ,&
            muip1jm1(grid%x_grd_num_u   ,grid%y_grd_num_u)  ,&
            muim1jp1(grid%x_grd_num_u   ,grid%y_grd_num_u)  ,&
            muim1jm1(grid%x_grd_num_u   ,grid%y_grd_num_u)  ,&
    
            vip1    (grid%x_grd_num_v   ,grid%y_grd_num_v)  ,&
            vim1    (grid%x_grd_num_v   ,grid%y_grd_num_v)  ,&
            vjp1    (grid%x_grd_num_v   ,grid%y_grd_num_v)  ,&
            vjm1    (grid%x_grd_num_v   ,grid%y_grd_num_v)  ,&
            vip1jp1 (grid%x_grd_num_v   ,grid%y_grd_num_v)  ,&
            vip1jm1 (grid%x_grd_num_v   ,grid%y_grd_num_v)  ,&
            vim1jp1 (grid%x_grd_num_v   ,grid%y_grd_num_v)  ,&
            vim1jm1 (grid%x_grd_num_v   ,grid%y_grd_num_v)  ,&
            mvip1   (grid%x_grd_num_v   ,grid%y_grd_num_v)  ,&
            mvim1   (grid%x_grd_num_v   ,grid%y_grd_num_v)  ,&
            mvjp1   (grid%x_grd_num_v   ,grid%y_grd_num_v)  ,&
            mvjm1   (grid%x_grd_num_v   ,grid%y_grd_num_v)  ,&
            mvip1jp1(grid%x_grd_num_v   ,grid%y_grd_num_v)  ,&
            mvip1jm1(grid%x_grd_num_v   ,grid%y_grd_num_v)  ,&
            mvim1jp1(grid%x_grd_num_v   ,grid%y_grd_num_v)  ,&
            mvim1jm1(grid%x_grd_num_v   ,grid%y_grd_num_v)  ,&
            
            zip1    (grid%x_grd_num     ,grid%y_grd_num)    ,&
            zim1    (grid%x_grd_num     ,grid%y_grd_num)    ,&
            zjp1    (grid%x_grd_num     ,grid%y_grd_num)    ,&
            zjm1    (grid%x_grd_num     ,grid%y_grd_num)    ,&
            zip1jp1 (grid%x_grd_num     ,grid%y_grd_num)    ,&
            zip1jm1 (grid%x_grd_num     ,grid%y_grd_num)    ,&
            zim1jp1 (grid%x_grd_num     ,grid%y_grd_num)    ,&
            zim1jm1 (grid%x_grd_num     ,grid%y_grd_num)    ,&
            mip1    (grid%x_grd_num     ,grid%y_grd_num)    ,&
            mim1    (grid%x_grd_num     ,grid%y_grd_num)    ,&
            mjp1    (grid%x_grd_num     ,grid%y_grd_num)    ,&
            mjm1    (grid%x_grd_num     ,grid%y_grd_num)    ,&
            mip1jp1 (grid%x_grd_num     ,grid%y_grd_num)    ,&
            mip1jm1 (grid%x_grd_num     ,grid%y_grd_num)    ,&
            mim1jp1 (grid%x_grd_num     ,grid%y_grd_num)    ,&
            mim1jm1 (grid%x_grd_num     ,grid%y_grd_num)
            
    real    zmip1   (grid%x_grd_num     ,grid%y_grd_num)    ,&
            zmim1   (grid%x_grd_num     ,grid%y_grd_num)    ,&
            zmjp1   (grid%x_grd_num     ,grid%y_grd_num)    ,&
            zmjm1   (grid%x_grd_num     ,grid%y_grd_num)    ,&
            zmip1jp1(grid%x_grd_num     ,grid%y_grd_num)    ,&
            zmip1jm1(grid%x_grd_num     ,grid%y_grd_num)    ,&
            zmim1jp1(grid%x_grd_num     ,grid%y_grd_num)    ,&
            zmim1jm1(grid%x_grd_num     ,grid%y_grd_num)
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
    
    !Matrix Variables
    !$omp parallel
    !$omp sections
    !$omp section
    uip1   (1:grid%x_grd_num_u-1,:)                     =   u   (2:grid%x_grd_num_u,:)      
    !$omp section
    uim1   (2:grid%x_grd_num_u,:)                       =   u   (1:grid%x_grd_num_u-1,:)     
    !$omp section
    ujp1   (:,1:grid%y_grd_num_u-1)                     =   u   (:,2:grid%y_grd_num_u)                     
    !$omp section
    ujm1   (:,2:grid%y_grd_num_u)                       =   u   (:,1:grid%y_grd_num_u-1)                   
    !!$omp section
    !uip1jp1(1:grid%x_grd_num_u-1,1:grid%y_grd_num_u-1)  =   u   (2:grid%x_grd_num_u  ,2:grid%y_grd_num_u)  
    !$omp section
    uip1jm1(1:grid%x_grd_num_u-1,2:grid%y_grd_num_u)    =   u   (2:grid%x_grd_num_u  ,1:grid%y_grd_num_u-1)
    !!$omp section
    !uim1jp1(2:grid%x_grd_num_u  ,1:grid%y_grd_num_u-1)  =   u   (1:grid%x_grd_num_u-1,2:grid%y_grd_num_u)  
    !!$omp section
    !uim1jm1(2:grid%x_grd_num_u  ,2:grid%y_grd_num_u)    =   u   (1:grid%x_grd_num_u-1,1:grid%y_grd_num_u-1)
    !$omp section
    muip1   (1:grid%x_grd_num_u-1,:)                    =   mu  (2:grid%x_grd_num_u,:)                     
    !$omp section
    muim1   (2:grid%x_grd_num_u,:)                      =   mu  (1:grid%x_grd_num_u-1,:)                   
    !!$omp section
    !mujp1   (:,1:grid%y_grd_num_u-1)                    =   mu  (:,2:grid%y_grd_num_u)                     
    !$omp section
    mujm1   (:,2:grid%y_grd_num_u)                      =   mu  (:,1:grid%y_grd_num_u-1)                   
    !!$omp section
    !muip1jp1(1:grid%x_grd_num_u-1,1:grid%y_grd_num_u-1) =   mu  (2:grid%x_grd_num_u  ,2:grid%y_grd_num_u)  
    !$omp section
    muip1jm1(1:grid%x_grd_num_u-1,2:grid%y_grd_num_u)   =   mu  (2:grid%x_grd_num_u  ,1:grid%y_grd_num_u-1)
    !!$omp section
    !muim1jp1(2:grid%x_grd_num_u  ,1:grid%y_grd_num_u-1) =   mu  (1:grid%x_grd_num_u-1,2:grid%y_grd_num_u)  
    !!$omp section
    !muim1jm1(2:grid%x_grd_num_u  ,2:grid%y_grd_num_u)   =   mu  (1:grid%x_grd_num_u-1,1:grid%y_grd_num_u-1)
    
    !$omp section
    vip1   (1:grid%x_grd_num_v-1,:)                     =   v   (2:grid%x_grd_num_v,:)                     
    !$omp section
    vim1   (2:grid%x_grd_num_v,:)                       =   v   (1:grid%x_grd_num_v-1,:)                   
    !$omp section
    vjp1   (:,1:grid%y_grd_num_v-1)                     =   v   (:,2:grid%y_grd_num_v)                     
    !$omp section
    vjm1   (:,2:grid%y_grd_num_v)                       =   v   (:,1:grid%y_grd_num_v-1)                   
    !!$omp section
    !vip1jp1(1:grid%x_grd_num_v-1,1:grid%y_grd_num_v-1)  =   v   (2:grid%x_grd_num_v  ,2:grid%y_grd_num_v)  
    !!$omp section
    !vip1jm1(1:grid%x_grd_num_v-1,2:grid%y_grd_num_v)    =   v   (2:grid%x_grd_num_v  ,1:grid%y_grd_num_v-1)
    !$omp section
    vim1jp1(2:grid%x_grd_num_v  ,1:grid%y_grd_num_v-1)  =   v   (1:grid%x_grd_num_v-1,2:grid%y_grd_num_v)
    !!$omp section
    !vim1jm1(2:grid%x_grd_num_v  ,2:grid%y_grd_num_v)    =   v   (1:grid%x_grd_num_v-1,1:grid%y_grd_num_v-1)
    !!$omp section
    !mvip1   (1:grid%x_grd_num_v-1,:)                    =   mv  (2:grid%x_grd_num_v,:)                     
    !!$omp section
    !mvim1   (2:grid%x_grd_num_v,:)                      =   mv  (1:grid%x_grd_num_v-1,:)                   
    !$omp section
    mvjp1   (:,1:grid%y_grd_num_v-1)                    =   mv  (:,2:grid%y_grd_num_v)                     
    !$omp section
    mvjm1   (:,2:grid%y_grd_num_v)                      =   mv  (:,1:grid%y_grd_num_v-1)                   
    !!$omp section
    !mvip1jp1(1:grid%x_grd_num_v-1,1:grid%y_grd_num_v-1) =   mv  (2:grid%x_grd_num_v  ,2:grid%y_grd_num_v)  
    !!$omp section
    !mvip1jm1(1:grid%x_grd_num_v-1,2:grid%y_grd_num_v)   =   mv  (2:grid%x_grd_num_v  ,1:grid%y_grd_num_v-1)
    !$omp section
    mvim1jp1(2:grid%x_grd_num_v  ,1:grid%y_grd_num_v-1) =   mv  (1:grid%x_grd_num_v-1,2:grid%y_grd_num_v)  
    !!$omp section
    !mvim1jm1(2:grid%x_grd_num_v  ,2:grid%y_grd_num_v)   =   mv  (1:grid%x_grd_num_v-1,1:grid%y_grd_num_v-1)
    
    !!$omp section
    !zip1    (1:grid%x_grd_num-1,:)                      =   z   (2:grid%x_grd_num,:)                    
    !$omp section
    zim1    (2:grid%x_grd_num,:)                        =   z   (1:grid%x_grd_num-1,:)                  
    !!$omp section
    !zjp1    (:,1:grid%y_grd_num-1)                      =   z   (:,2:grid%y_grd_num)                    
    !$omp section
    zjm1    (:,2:grid%y_grd_num)                        =   z   (:,1:grid%y_grd_num-1)                  
    !!$omp section
    !zip1jp1 (1:grid%x_grd_num-1,1:grid%y_grd_num-1)     =   z   (2:grid%x_grd_num  ,2:grid%y_grd_num)   
    !!$omp section
    !zip1jm1 (1:grid%x_grd_num-1,2:grid%y_grd_num)       =   z   (2:grid%x_grd_num  ,1:grid%y_grd_num-1) 
    !!$omp section
    !zim1jp1 (2:grid%x_grd_num  ,1:grid%y_grd_num-1)     =   z   (1:grid%x_grd_num-1,2:grid%y_grd_num)   
    !!$omp section
    !zim1jm1 (2:grid%x_grd_num  ,2:grid%y_grd_num)       =   z   (1:grid%x_grd_num-1,1:grid%y_grd_num-1)
    !!$omp section
    !mip1   (1:grid%x_grd_num-1,:)                       =   m   (2:grid%x_grd_num,:)                    
    !!$omp section
    !mim1   (2:grid%x_grd_num,:)                         =   m   (1:grid%x_grd_num-1,:)                  
    !!$omp section
    !mjp1   (:,1:grid%y_grd_num-1)                       =   m   (:,2:grid%y_grd_num)                    
    !!$omp section
    !mjm1   (:,2:grid%y_grd_num)                         =   m   (:,1:grid%y_grd_num-1)                  
    !!$omp section
    !mip1jp1(1:grid%x_grd_num-1,1:grid%y_grd_num-1)      =   m   (2:grid%x_grd_num  ,2:grid%y_grd_num)   
    !!$omp section
    !mip1jm1(1:grid%x_grd_num-1,2:grid%y_grd_num)        =   m   (2:grid%x_grd_num  ,1:grid%y_grd_num-1) 
    !!$omp section
    !mim1jp1(2:grid%x_grd_num  ,1:grid%y_grd_num-1)      =   m   (1:grid%x_grd_num-1,2:grid%y_grd_num)   
    !!$omp section
    !mim1jm1(2:grid%x_grd_num  ,2:grid%y_grd_num)        =   m   (1:grid%x_grd_num-1,1:grid%y_grd_num-1)
    
    !$omp section
    zmip1    (1:grid%x_grd_num-1,:)                     =   zm  (2:grid%x_grd_num,:)                    
    !$omp section
    zmim1    (2:grid%x_grd_num,:)                       =   zm  (1:grid%x_grd_num-1,:)                  
    !$omp section
    zmjp1    (:,1:grid%y_grd_num-1)                     =   zm  (:,2:grid%y_grd_num)                    
    !$omp section
    zmjm1    (:,2:grid%y_grd_num)                       =   zm  (:,1:grid%y_grd_num-1)                  
    !!$omp section
    !zmip1jp1 (1:grid%x_grd_num-1,1:grid%y_grd_num-1)    =   zm  (2:grid%x_grd_num  ,2:grid%y_grd_num)   
    !!$omp section
    !zmip1jm1 (1:grid%x_grd_num-1,2:grid%y_grd_num)      =   zm  (2:grid%x_grd_num  ,1:grid%y_grd_num-1) 
    !!$omp section
    !zmim1jp1 (2:grid%x_grd_num  ,1:grid%y_grd_num-1)    =   zm  (1:grid%x_grd_num-1,2:grid%y_grd_num)   
    !!$omp section
    !zmim1jm1 (2:grid%x_grd_num  ,2:grid%y_grd_num)      =   zm  (1:grid%x_grd_num-1,1:grid%y_grd_num-1)
    
    !$omp end sections
    !$omp end parallel
    
    !Matrix solve
    !$omp parallel
    !$omp sections
    !$omp section
    u_t=-mu*c*((uip1+u)      *(uip1-u)        &   
              +(u+uim1)      *(u-uim1)        &   
              +(vjp1+vim1jp1)*(ujp1-u)        &
              +(v+vim1)      *(u-ujm1)        &   
              +4.d0*g        *(z-zim1))       &            
        +0.25*(fu                             &     
              +c*((uip1+u)      *(mvjp1-mv)   &      
                 +(u+uim1)      *(mvim1jp1-mv)&       
                 -(vjp1+v)      *(muip1-mu)   &      
                 -(vim1jp1+vim1)*(mu-muim1))) &  
             *(vjp1+vim1jp1+v+vim1)
    u_t(:,1)=0
    u_t(1,:)=0
    u_t(:,grid%y_grd_num_u)=0
    u_t(grid%x_grd_num_u,:)=0
    
    !$omp section
    v_t=-mv*c*((uip1+uip1jm1)*(vip1-v)       &
              +(u+ujm1)*(v-vim1)             &
              +(vjp1+v)*(vjp1-v)             &
              +(v+vjm1)*(v-vjm1)             &
              +4.d0*g*(z-zjm1))              &   
        -0.25*(fv                            & 
              +c*((uip1+u)*(mvjp1-mv)        &
                 +(uip1jm1+ujm1)*(mv-mvjm1)  &
                 -(vjp1+v)*(muip1-mu)        &
                 -(v+vjm1)*(muip1jm1-mujm1)))&
             *(uip1+u+uip1jm1+ujm1)
    v_t(:,1)=0
    v_t(1,:)=0
    v_t(:,grid%y_grd_num_v)=0
    v_t(grid%x_grd_num_v,:)=0
    
    !$omp section
    z_t=-m**2*2.d0*c*(uip1*(zmip1-zm)& 
                     +u   *(zm-zmim1)&        
                     +vjp1*(zmjp1-zm)&      
                     +v   *(zm-zmjm1)&        
                     +2.d0*(zm-z0)   &              
                    *(uip1-u+vjp1-v))
    z_t(:,1)=0
    z_t(1,:)=0
    z_t(:,grid%y_grd_num)=0
    z_t(grid%x_grd_num,:)=0
    !$omp end sections
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

    end subroutine time_integration_matrix

    
    
    subroutine time_integration_loop(ua,va,za,u,v,z,uc,vc,zc,ts,option)
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

    end subroutine time_integration_loop
    
end module integration