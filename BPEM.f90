    !Barotropic primitive equation model
    !Written by Zhou Lilong Dec. 2017

    program Barotropic_Model
    use control
    use subroutines
    use initialization
    use solve
    use output
    implicit none
    integer t1,t2

    call SYSTEM_CLOCK(t1)

    call read_namelist
    call choose_slash

    !Input
    call init_data

    !Open output file
    call create_output

    call Runge_Kutta3

    call SYSTEM_CLOCK(t2)
    print*,'It took ',dble(t2-t1)/time_divide_value,' seconds to run this program'
    end program Barotropic_Model
