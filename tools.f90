module tools
    contains
    !���ӳ�����չ�������գ��ֿ��Լ���������ʱ����仯���������ʱ����
    !character*16 option
    !integer value
    !option='calculate_second'�������
    !option='calculate_minute'���ּ���
    !option='calculate_hour'��ʱ����
    !option='calculate_day'���ռ���
    !valueΪ��������
    !��Ҫ��10������value=10��option='calculate_minute'��������value=-8��option='calculate_second'
    subroutine time_calculator(year_in,month_in,day_in,hour_in,minute_in,second_in,&
                               year_out,month_out,day_out,hour_out,minute_out,second_out,&
                               option,value)
        implicit none
        integer,intent(in)::year_in,month_in,day_in,hour_in,minute_in,second_in
        integer,intent(out),optional::year_out,month_out,day_out,hour_out,minute_out,second_out
        integer,intent(in)::value
        character*(*),intent(in)::option

        integer year_temp,month_temp,day_temp,hour_temp,minute_temp,second_temp
        integer day_increase,hour_increase,minute_increase!,second_increase

        !Julian variables
        real*8 Julian_day
        real*8 temp1,temp2,temp3
        integer::J,N,L1,Y1,L2,M1,L3

        year_temp   = year_in
        month_temp  = month_in
        day_temp    = day_in
        hour_temp   = hour_in
        minute_temp = minute_in
        second_temp = second_in

        minute_increase=0
        hour_increase=0
        day_increase=0

        select case(trim(adjustl(option)))
        case('calculate_second')
            second_temp=second_in+value
        case('calculate_minute')
            minute_temp=minute_in+value
        case('calculate_hour')
            hour_temp=hour_in+value
        case('calculate_day')
            day_temp=day_in+value
        end select

        !��仯
        if(second_temp>=0)then
            minute_increase=minute_increase+int(second_temp/60)
            if(present(second_out))second_out=mod(second_temp,60)
        else
            minute_increase=minute_increase-1+int(second_temp/60)
            if(present(second_out))second_out=60+mod(second_temp,60)
        endif

        !���ӱ仯
        minute_temp=minute_temp+minute_increase
        if(minute_temp>=0)then
            hour_increase=hour_increase+int(minute_temp/60)
            if(present(minute_out))minute_out=mod(minute_temp,60)
        else
            hour_increase=hour_increase-1+int(minute_temp/60)
            if(present(minute_out))minute_out=60+mod(minute_temp,60)
        endif

        !Сʱ�仯
        hour_temp=hour_temp+hour_increase
        if(hour_temp>=0)then
            day_increase=day_increase+int(hour_temp/24)
            if(present(hour_out))hour_out=mod(hour_temp,24)
        else
            day_increase=day_increase-1+int(hour_temp/24)
            if(present(hour_out))hour_out=24+mod(hour_temp,24)
        endif

        !�����ձ仯
        !����ʱ��仯���������
        temp1=1461*(year_in+4800+(month_in-14)/12)/4
        temp2=367*(month_in-2-((month_in-14)/12)*12)/12
        temp3=3*(year_in+4900+(month_in-14)/12)/100/4

        Julian_day=day_in-32075+temp1+temp2-temp3-0.5d0+day_increase

        !����ʱ��仯�����ʵ������
        J=int(Julian_day+0.5D0)
        N=4*(J+68569)/146097
        L1=J+68569-(N*146097+3)/4
        Y1=4000*(L1+1)/1461001
        L2=L1-1461*Y1/4+31

        M1=80*L2/2447
        if(present(day_out))day_out=L2-2447*M1/80
        L3=INT(M1/11)
        if(present(month_out))month_out=M1+2-12*L3
        if(present(year_out))year_out=100*(N-49)+Y1+L3

    end subroutine time_calculator

    !���Ӻ���Ϊnetcdf�Ĵ����������Է��ش���nc�ļ�ʱ�Ĵ�����Ϣ
    subroutine handle_err(status)
        use netcdf
        !use netcdf90
        implicit none
        include 'netcdf.inc'
        integer, intent (in) :: status

        if(status /= nf_noerr) then
            print *, trim(nf_strerror(status))
            stop "Stopped"
        end if
    end subroutine handle_err

    !���Ӻ����ɽ����α���ת��Ϊ�ַ���
    function num2str(intNum)
        integer,intent(in) :: intNum
        character(len=:),allocatable:: num2str
        character(len=255):: tmpstr
        write(tmpstr,*)intNum
        allocate(character(len=len_trim(adjustl(tmpstr)))::num2str)
        num2str=trim(adjustl(tmpstr))
    end function num2str

    !���Ӻ����ɽ��ַ���ת��Ϊ˫�������ͱ���
    function str2num(CharNum)
        character*(*),intent(in) :: CharNum
        real*8,allocatable:: str2num
        read(CharNum,*) str2num
    end function str2num

    !���ӳ������ڶ�ȡMicaps��diamond3��ʽ����
    !������վ��(integer,id(station_number))����γ��(real lon(station_number),lat(station_number))������(real altitude(station_number))����ˮ��(real rain(station_number))
    !��Ҫ����file�����ļ�·���Լ��ļ�����skip�ӵ�һ�п�ʼ������������վ������������������Ҫ������������һ��Ϊ11�У���������һ�м�¼�����б�Ե���ϸ���ľ�γ�ȼ�¼����վ������
    subroutine M3_read(file,unit,skip,station_number,id,lon,lat,altitude,rain)
        implicit none
        character*200,intent(in) :: file
        integer,intent(in) :: unit,skip
        integer,intent(inout) ::station_number
        integer,intent(out),optional,allocatable :: id(:),altitude(:)
        real,intent(out),optional,allocatable :: lon(:),lat(:),rain(:)

        integer i,sig,iostat,count_station_number

        open(unit,file=trim(file),status='old')
        !���������ļ�ͷ
        do i=1,skip
            read(unit,*)
        enddo

        if(station_number<=0)then
            iostat=0
            count_station_number=0
            do while(iostat==0)
                read(unit,*,iostat=iostat)
                if(iostat==0)then
                    count_station_number=count_station_number+1
                else
                    print*,'�ļ� ',trim(adjustl(file)),' ��վ�����Ϊ��',count_station_number
                    station_number=count_station_number
                endif
            enddo
        endif

        if(present(id).or.present(lon).or.present(lat).or.present(altitude).or.present(rain))then
            allocate(id(station_number),lon(station_number),lat(station_number),altitude(station_number),rain(station_number))
            count_station_number=0
            do i=1,station_number
                read(unit,*,iostat=iostat)id(i),lon(i),lat(i),altitude(i),rain(i)
            enddo
        endif

        close(unit)
    end subroutine M3_read


    subroutine M8_read(file,unit,station_number,id,lon,lat,altitude,&
                       weather1,wind_direction1,wind_speed1,&
                       minimum_temperature,maximum_temperature,&
                       weather2,wind_direction2,wind_speed2)
        implicit none
        character*(*),intent(in) :: file
        integer,intent(in) :: unit
        integer,intent(inout) ::station_number
        integer,intent(out),optional,allocatable :: id(:)
        real,intent(out),optional,allocatable :: lon(:),lat(:)
        integer,intent(out),optional,allocatable:: altitude(:),weather1(:),wind_direction1(:),wind_speed1(:),&
                                                   weather2(:),wind_direction2(:),wind_speed2(:)
        real,intent(out),optional,allocatable :: minimum_temperature(:),maximum_temperature(:)
        integer i,iostat,count_station_number

        open(unit,file=trim(adjustl(file)),status='old')
        read(unit,*)
        read(unit,*)
        if(station_number<=0)then
            iostat=0
            count_station_number=0
            do while(iostat==0)
                read(unit,*,iostat=iostat)
                if(iostat==0)then
                    count_station_number=count_station_number+1
                else
                    print*,'�ļ� ',trim(adjustl(file)),' ��վ�����Ϊ��',count_station_number
                    station_number=count_station_number
                endif
            enddo
            rewind(unit)
            read(unit,*)
        else
            read(unit,*)
        endif

        if(present(id)                  .or.present(lon)                .or.present(lat)        .or.present(altitude).or.&
           present(weather1)            .or.present(wind_direction1)    .or.present(wind_speed1).or.&
           present(minimum_temperature) .or.present(maximum_temperature).or.&
           present(weather2)            .or.present(wind_direction2)    .or.present(wind_speed2))then
            allocate(id(station_number),&
                     lon(station_number),&
                     lat(station_number),&
                     altitude(station_number),&
                     weather1(station_number),&
                     wind_direction1(station_number),&
                     wind_speed1(station_number),&
                     minimum_temperature(station_number),&
                     maximum_temperature(station_number),&
                     weather2(station_number),&
                     wind_direction2(station_number),&
                     wind_speed2(station_number))

            do i=1,station_number
                read(unit,*,iostat=iostat)id(i),lon(i),lat(i),altitude(i),&
                                          weather1(i),wind_direction1(i),wind_speed1(i),&
                                          minimum_temperature(i),maximum_temperature(i),&
                                          weather2(i),wind_direction2(i),wind_speed2(i)
            enddo
        endif
        close(unit)
    end subroutine M8_read

    !дMicaps��3���ļ�
    !level=-1 ��ʾ��6Сʱ��ˮ��������ˮ��Ϊ0.0mmʱ��T������ˮ��Ϊ0.1-0.9ʱ��һλС��������ˮ������1ʱֻ��������
    !level=-2 ��ʾ��24Сʱ��ˮ��������ˮ��С��1mmʱ������ڵ���1mmʱֻ��������
    !level=-3 ��ʾ���¶ȣ�ֻ������
    !level��������ֵ����ԭ����
    subroutine M3_write(file,unit,data_title,station_number,station_id,longitude,latitude,altitude,rainfall,level,undefined_value,year,month,day,hour)
        implicit none
        character*(*),intent(in):: file
        integer,intent(in):: unit,station_number
        integer,intent(in):: station_id(station_number)
        real,intent(in):: longitude(station_number),latitude(station_number),rainfall(station_number)
        integer,intent(in),optional:: altitude(station_number)
        character*(*),intent(in),optional:: data_title
        integer,intent(in),optional:: level
        real,intent(in):: undefined_value
        character*2,intent(in),optional:: year,month,day,hour

        integer i
        integer sta_n
        real alt(station_number)
        character*2 yy,mm,dd,hh
        character*200 title
        integer lev

        if(.not.present(data_title))then
            title='δ֪����'
        else
            title=data_title
        endif

        if((.not.present(year)).or.(.not.present(month)).or.(.not.present(day)).or.(.not.present(hour)))then
            yy='17'
            mm='01'
            dd='01'
            hh='00'
        else
            yy=year
            mm=month
            dd=day
            hh=hour
        endif

        if(.not.present(level))then
            lev=1000
        else
            lev=level
        endif

        if(.not.present(altitude))then
            alt=9999
        else
            alt=altitude
        endif

        sta_n=count(longitude>=-180.and.longitude<=180.and.latitude>=-90.and.latitude<=90.and.abs(altitude)<=10000.and.rainfall/=undefined_value)

        open(unit,file=trim(adjustl(file)))
        write(unit,'(a)')'diamond 3 '//title
        write(unit,'(4(a2,1x),i6,a12)')yy,mm,dd,hh,lev,'  0  0  0  0'
        write(unit,'(a2,i10)')'1 ',sta_n
        do i=1,station_number
            if(longitude(i)>=-180.and.longitude(i)<=180.and.latitude(i)>=-90.and.latitude(i)<=90.and.abs(altitude(i))<=10000.and.rainfall(i)/=undefined_value)write(unit,'(i7,1x,f6.2,1x,f6.2,1x,f7.1,1x,f8.2)')station_id(i),longitude(i),latitude(i),alt(i),rainfall(i)
        enddo
        close(unit)
    end subroutine M3_write

    logical function usable_judgement(array,element_num,upper_bound,lower_bound)
        implicit none
        integer,intent(in):: element_num
        real,intent(in):: array(element_num)
        real,intent(in):: upper_bound,lower_bound
        integer i
        integer undef_num
        real forward_gradient(element_num),&    !ǰ���ݶ�
             backward_gradient(element_num),&   !����ݶ�
             central_gradient(element_num)      !������ݶ�

        usable_judgement=.true.

        if(element_num>2)then
            !���ȱ��ֵռ�ȳ���1/4��ֱ���ж�Ϊ��Ч������
            undef_num=count(array<lower_bound.and.array>upper_bound)
            if(undef_num/element_num>=0.25)then
                usable_judgement=.false.
                return
            endif

            !ǰ���ʽ�����ݶ�
            do i=1,element_num-1
                if(array(i)>=lower_bound.and.array(i)<=upper_bound.and.&
                   array(i+1)>=lower_bound.and.array(i+1)<=upper_bound)then
                    forward_gradient(i)=array(i+1)-array(i)
                else
                    forward_gradient(i)=9999
                endif
            enddo

            !����ʽ�����ݶ�
            do i=2,element_num
                if(array(i-1)>=lower_bound.and.array(i-1)<=upper_bound.and.&
                   array(i)>=lower_bound.and.array(i)<=upper_bound)then
                    backward_gradient(i)=array(i)-array(i-1)
                else
                    backward_gradient(i)=9999
                endif
            enddo

            !������ʽ�����ݶ�
            do i=2,element_num-1
                if(array(i-1)>=lower_bound.and.array(i-1)<=upper_bound.and.&
                   array(i+1)>=lower_bound.and.array(i+1)<=upper_bound)then
                    central_gradient(i)=(array(i+1)-array(i-1))/2
                else
                    central_gradient(i)=9999
                endif
            enddo

            do i=2,element_num-1
                if(central_gradient(i)/=9999.and.&
                   forward_gradient(i)==9999.and.&
                   backward_gradient(i)==9999)then
                    if(backward_gradient(i-1)/=9999.and.forward_gradient(i+1)/=9999)then
                        if(backward_gradient(i-1)/abs(backward_gradient(i-1))/=forward_gradient(i+1)/abs(forward_gradient(i+1)))then
                            usable_judgement=.false.
                            return
                        endif
                    else
                       usable_judgement=.false.
                       return
                    endif
                endif
            enddo
        endif
    end function usable_judgement

end module tools