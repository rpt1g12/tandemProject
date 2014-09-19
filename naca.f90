 program naca
    implicit none
    integer,parameter :: nr=kind(0.0d0)
    integer,parameter :: ntot=1000
    real(nr),parameter :: pi=3.141592653589793_nr
    integer :: i
    real(kind=nr) :: c,t,k
    real(kind=nr) :: r,h
    real(kind=nr), dimension(0:ntot)  ::  y,x,xc

    c = 1; t = 12;
    t=t*0.01e0
    h = (pi/real(ntot))
    k=0.991148635e0
    open (unit=1, file='aerofoil.dat')

    do i = 0, ntot
       r = i*h
       xc(i) = 0.5_nr*(1_nr-cos(r))
       y(i) =(t*c*k/0.2e0) &
       *(0.298222773e0*sqrt(xc(i))-0.127125232e0*xc(i)- 0.357907906e0*xc(i)**2+0.291984971e0*xc(i)**3-0.105174606e0*xc(i)**4)
    if (i==ntot) then
       y(i)=0e0
    end if
    write(1,'(f9.7,2es15.7)') xc(i),-y(i)
    end do

    do i = 0, ntot
    write(1,'(f9.7,2es15.7)') xc(i),y(i)
    end do

    close(1)
 end program naca
