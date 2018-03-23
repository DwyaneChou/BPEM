    module solve
    use control
    use subroutines
    use integration
    use output
    use omp_lib
    contains
    !3rd order Runge Kutta time integeration
    subroutine Runge_Kutta3
    implicit none
    type runge_kutta
        real,allocatable:: u(:,:),v(:,:),z(:,:)
    endtype
    type(runge_kutta)::rk1,rk2,rk3
    real ua_temp(grid%x_grd_num_u,grid%y_grd_num_u),va_temp(grid%x_grd_num_v,grid%y_grd_num_v),za_temp(grid%x_grd_num,grid%y_grd_num),&
         ub_temp(grid%x_grd_num_u,grid%y_grd_num_u),vb_temp(grid%x_grd_num_v,grid%y_grd_num_v),zb_temp(grid%x_grd_num,grid%y_grd_num),&
         uc_temp(grid%x_grd_num_u,grid%y_grd_num_u),vc_temp(grid%x_grd_num_v,grid%y_grd_num_v),zc_temp(grid%x_grd_num,grid%y_grd_num)

    integer i,loop_num
    integer::output_num=1
    real integrated_step_num

    allocate(rk1%u(grid%x_grd_num_u,grid%y_grd_num_u),rk1%v(grid%x_grd_num_v,grid%y_grd_num_v),rk1%z(grid%x_grd_num,grid%y_grd_num),&
             rk2%u(grid%x_grd_num_u,grid%y_grd_num_u),rk2%v(grid%x_grd_num_v,grid%y_grd_num_v),rk2%z(grid%x_grd_num,grid%y_grd_num),&
             rk3%u(grid%x_grd_num_u,grid%y_grd_num_u),rk3%v(grid%x_grd_num_v,grid%y_grd_num_v),rk3%z(grid%x_grd_num,grid%y_grd_num))

    !Calculate Lateral Boundary Condition Weight Function Coefficient
    call LBC_WFC(grid%F1,grid%F2,grid%bdy_width,grid%Spec_Zone,grid%Relax_Zone,grid%dx)
    
    loop_num=ceiling(run_hours*3600/time_step)
    
    grid%ub = grid%ua
    grid%vb = grid%va
    grid%zb = grid%za
    
    do i=1,loop_num
        grid%integerated_step_num=i

        call diffusion(grid%ub,grid%x_grd_num_u ,grid%y_grd_num_u)
        call diffusion(grid%vb,grid%x_grd_num_v ,grid%y_grd_num_v)
        call diffusion(grid%zb,grid%x_grd_num   ,grid%y_grd_num)
        
        call RK3_WRF(grid%ub,grid%vb,grid%zb,grid%uc,grid%vc,grid%zc,time_step)

        call ULBC_WRF(grid%uc,grid%vc,grid%zc,grid%integerated_step_num,time_step)
        
        !Output data
        if(mod(time_step*i,history_interval_seconds)==0.or.i==loop_num)then
            output_num=output_num+1
            print*,'It is time to write the output at ',grid%times(:,output_num)
            call write_output(output_num,grid%za,grid%ua,grid%va)
        endif
        
        ! Transmit arrays for the next loop
        call transmit_arrays(grid%ua,grid%va,grid%za,grid%ub,grid%vb,grid%zb)
        call transmit_arrays(grid%ub,grid%vb,grid%zb,grid%uc,grid%vc,grid%zc)
        
    enddo

    end subroutine Runge_Kutta3

    !RK3 from WRF
    subroutine RK3_WRF(u_in,v_in,z_in,u_out,v_out,z_out,ts)
    implicit none
    real,intent(in)::u_in(grid%x_grd_num_u,grid%y_grd_num_u),v_in(grid%x_grd_num_v,grid%y_grd_num_v),z_in(grid%x_grd_num,grid%y_grd_num)
    real,intent(in)::ts
    real,intent(out)::u_out(grid%x_grd_num_u,grid%y_grd_num_u),v_out(grid%x_grd_num_v,grid%y_grd_num_v),z_out(grid%x_grd_num,grid%y_grd_num)

    real ua(grid%x_grd_num_u,grid%y_grd_num_u),va(grid%x_grd_num_v,grid%y_grd_num_v),za(grid%x_grd_num,grid%y_grd_num),&
         ub(grid%x_grd_num_u,grid%y_grd_num_u),vb(grid%x_grd_num_v,grid%y_grd_num_v),zb(grid%x_grd_num,grid%y_grd_num),&
         uc(grid%x_grd_num_u,grid%y_grd_num_u),vc(grid%x_grd_num_v,grid%y_grd_num_v),zc(grid%x_grd_num,grid%y_grd_num),&
         ut(grid%x_grd_num_u,grid%y_grd_num_u),vt(grid%x_grd_num_v,grid%y_grd_num_v),zt(grid%x_grd_num,grid%y_grd_num)
    real dt_times,sub_time_step1,sub_time_step2,sub_time_step3
    integer loop_num
    integer i

    dt_times=int(ts)/int(time_step)
    if(dt_times<=1)then
        sub_time_step1=dble(ts)
        loop_num=1
    elseif(dt_times>1)then
        sub_time_step1=ts/dble(loop_num)
        loop_num=ceiling(dt_times)
    endif

    sub_time_step3=sub_time_step1/3.d0
    sub_time_step2=sub_time_step1/2.d0
    
    ua=u_in
    va=v_in
    za=z_in
    
    do i=1,loop_num
        call time_integration(ua,va,za,ua,va,za,ut,vt,zt,sub_time_step3,2)
        ub=ua+ut
        vb=va+vt
        zb=za+zt
        
        call time_integration(ua,va,za,ub,vb,zb,ut,vt,zt,sub_time_step2,2)
        uc=ua+ut
        vc=va+vt
        zc=za+zt
        
        call time_integration(ub,vb,zb,uc,vc,zc,ut,vt,zt,sub_time_step1,2)
        ua=ua+ut
        va=va+vt
        za=za+zt
    enddo

    u_out=ua
    v_out=va
    z_out=za
    end subroutine RK3_WRF
    
    end module solve