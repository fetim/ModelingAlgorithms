Program Modeling
!$ use omp_lib    
implicit none

integer :: k,i,j
integer :: count_snap,Nsnap
integer :: Nx,Nz,Nt,Nchannel,Nsamples, Nshot
integer :: shot, reg_snapshot
integer :: Nt_src
integer :: rx, rz
integer, allocatable,dimension(:) :: sx,sz

integer :: inix,iniz,endx,endz
real :: dx,dz,dt
real :: fcut, fcut_aux
real :: vpmax,vpmin
real :: p_xx, p_zz
real :: t_src,t0_src,src_aux
real, parameter :: pi = 4.0*ATAN(1.0)

real,allocatable,dimension(:) :: source, parameters
real,allocatable,dimension(:) :: P1,P2,P3,VP
real,allocatable,dimension(:) :: Seism, record
integer :: index

call importascii("../../parameters/2D_acoustic_modeling.dat",8,parameters)

! model parameters
Nx = int(parameters(1))
Nz = int(parameters(2))
dx = parameters(3)
dz = parameters(4)

write(*,*)"Nx = ", Nx
write(*,*)"Nz = ", Nz
write(*,*)"dx = ", dx
write(*,*)"dz = ", dz

Nchannel = Nx

inix=3
iniz=3
endx=Nx-1
endz=Nz-1

! time parameters
Nt = int(parameters(5))
dt = parameters(6)

write(*,*)"Nt = ", Nt
write(*,*)"dt = ", dt

Nsamples = Nt

!source parameters
Nshot = int(parameters(7))
fcut = parameters(8)

write(*,*)"Nshot = ", Nshot
write(*,*)"fcut = ", fcut

reg_snapshot = 1

! receiver depth
rz=20+1

! Allocate arrays
allocate(sx(Nshot))
allocate(sz(Nshot))
allocate(source(Nt))
allocate(VP(Nx*Nz))
allocate(P1(Nx*Nz))
allocate(P2(Nx*Nz))
allocate(P3(Nx*Nz))
allocate(Seism(Nsamples*Nchannel*Nshot))
allocate(record(Nsamples))

!Initializate arrays and counters
count_snap=1
Nsnap = 10
source = 0.0
P1 =0.0
P2 =0.0
P3 =0.0
Seism = 0.0
record = 0.0

do shot=1,Nshot
    sx(shot)   = Nx/2  + 5*(shot-1) + 1
    sz(shot)   = 20 + 1
end do

! velocity model
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Nx,Nz,VP)
!$OMP DO
    do index=1,Nx*Nz
        if (mod(index,Nz)==0) then
            i = index/Nz
        else    
            i = index/Nz + 1
        end if

        j = index  - (i-1)*Nz

        if (j<=100+1) then
            VP(index) = 1600
        else
            VP(index) = 1800
        end if

    end do
!$OMP END DO
!$OMP END PARALLEL

vpmax = maxval(VP)
vpmin = minval(VP)

call check_acoustic_stability(dt, dz, vpmax, vpmin, fcut, 2, 4) 

! Create Source Ricker
fcut_aux = fcut/(3.*sqrt(pi))          ! Ajust to cut of gaussian function
t0_src   = 4*sqrt(pi)/fcut             ! Initial time source
Nt_src   = nint(2*t0_src/dt) + 1       ! Number of elements of the source
do k=1,Nt_src                          ! Nts=nint(tm/dt)+1
    t_src=(k-1)*dt-t0_src                    !Delay Time
    src_aux=pi*(pi*fcut_aux*t_src)*(pi*fcut_aux*t_src)
    source(k) = -(2*src_aux-1)*exp(-src_aux)    
end do

! Register in disk
open(23, file='snapshots.bin', status='replace',&
&FORM='unformatted',ACCESS='direct', recl=(Nx*Nz*4))

open(24, file='seismogram.bin', status='replace',&
&FORM='unformatted',ACCESS='direct', recl=(Nsamples*Nchannel*Nshot*4))

do shot = 1,Nshot
    print*, "Running shot ...",shot
    ! Solve wave equation
    do k=1,Nt

        P2(sz(shot) + Nz*(sx(shot)-1)) = P2(sz(shot) + Nz*(sx(shot)-1)) - source(k)              
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(P1,P2,P3,VP,record,dt,dx,dz,Nx,Nz,inix,iniz,endx,endz,Nt,rz,Nchannel,Nsamples)
        !$OMP DO
        ! Vectorized spatial loop    
        do index=1,Nx*Nz 
            !4th order in space and 2nd order in time
            if (mod(index,Nz)==0) then
                i = index/Nz
            else    
                i = index/Nz + 1
            end if

            j = index  - (i-1)*Nz

            if ((i>=inix .and. j>=iniz) .and. (i<=endx .and. j<=endz)) then            
                p_zz = (-1*P2(j-2 + Nz*(i-1)) + 16*P2(j-1 + Nz*(i-1)) &
                - 30*P2(j + Nz*(i-1)) + 16*P2(j+1 + Nz*(i-1)) - 1*P2(j+2 + Nz*(i-1) ) )/(12*dz*dz)

                p_xx = (-1*P2(j + Nz*(i-1-2)) + 16*P2(j + Nz*(i-1-1)) &
                - 30*P2(j + Nz*(i-1)) + 16*P2(j + Nz*(i-1+1)) - 1*P2(j + Nz*(i-1+2) ) ) / (12*dx*dx)

                P3(index) = 2*P2(j + Nz*(i-1)) - P1(j + Nz*(i-1)) &
                + (dt*dt)*(VP(j + Nz*(i-1))*VP(j + Nz*(i-1)))*(p_xx + p_zz)

            end if
        end do            
        !$OMP END DO NOWAIT

        !Storage Seismogram        
        !$OMP DO
        do rx=1,Nchannel                               
            record(rx) = P3(rz + Nz*(rx-1))                
        end do
        !$OMP END DO NOWAIT      

        ! update fields                                
        !$OMP DO
        do index=1,Nx*Nz
            P1(index)=P2(index)
            P2(index)=P3(index)    
        end do
        !$OMP END DO NOWAIT       

        !$OMP END PARALLEL 

        !Storage Seismogram        
        do rx=1,Nchannel                               
            Seism(k + (rx-1)*Nsamples + (shot-1)*(Nchannel*Nsamples) ) = record(rx)
        end do

        !Register snapshots
        if ((mod((k-1),Nt/Nsnap) == 0) .and. (reg_snapshot == shot) ) then
            print*, "Propagation time =", (k-1)*dt, "Registering snapshot", count_snap
            write(23,rec=count_snap) ((P3(j + Nz*(i-1)) + 1.0e-3*VP(j + Nz*(i-1)),j=1,Nz),i=1,Nx)
            count_snap=count_snap+1
        end if     

    end do

    ! Restarting fields
    P1 =0.0
    P2 =0.0
    P3 =0.0

end do
!Register Seismogram
write(24,rec=1) (Seism(k),k=1,Nsamples*Nchannel*Nshot)

!close files
close(23)
close(24)

contains

SUBROUTINE importascii(name,N_lines,output)    
    IMPLICIT NONE
    LOGICAL                                        :: checkfile
    INTEGER                                        :: k
    CHARACTER(len=*),INTENT(in)                    :: name
    INTEGER,INTENT(in)                             :: N_lines
    REAL,DIMENSION(:),ALLOCATABLE,INTENT(out)      :: output

    write(*,*) "Reading ", name, " Number of lines = ", N_lines
    INQUIRE(file=name, exist=checkfile) !verify if wavelet file exist
    if (checkfile) then            
        ALLOCATE(output(N_lines))        
        open(78, file=name, status='unknown', form='formatted')               
        do k=1,N_lines
            read(78,*)output(k)
        end do
        close(78)
    else
        print*, 'Error opening file'           
        stop
    end if

    RETURN
END SUBROUTINE importascii

SUBROUTINE check_acoustic_stability(dt, dh, vpmax, vpmin, fcut, T_order, S_order)
    
    INTEGER, INTENT(IN) :: T_order, S_order
    REAL,    INTENT(IN) :: dt,dh,vpmax,vpmin,fcut

    REAL                :: alpha, beta, dh_max,dt_max


    if ((T_order == 2) .and. (S_order==4)) then	
        alpha  = 5
        beta   = 4
        dh_max = vpmin/(alpha*fcut)
        dt_max = dh_max/(beta*vpmax)

        if ( (dt <= dt_max) .and. (dh <= dh_max )) then

            print*, "Dispersion and stability conditions ... OK!"

        else
            print*,"Read time and space parameters"
            print*,"Time                        = ", dt
            print*,"Space                       = ", dh			

            print*,"Recommended time spacing and time interval"
            print*,"Time  - Stability limit  - dt < ", dh/(beta*vpmax)
            print*,"Time  - Stability limit  - dt < ", dt_max
            print*,"Space - Dispersion limit - dh < ", dh_max
            print*,"The program will be stopped!  "
            stop
        end if
    else	
        print*,"Select the right order "
        print*," Time derivative order = Space derivative order = , ", T_order, S_order
        stop
    end if

END SUBROUTINE check_acoustic_stability

end program