Program Modeling
!$ use omp_lib    
implicit none

integer :: k,i,j
integer :: count_snap, Nsnap
integer :: Nx,Nz,Nt
integer :: sx,sz,rx,rz
integer :: Nt_src

integer :: inix,iniz,endx,endz
real :: h,dt
real :: fcut, fcut_aux
real :: p_xx, p_zz
real :: t_src,t0_src,src_aux
real, parameter :: pi = 4.0*ATAN(1.0)

real,allocatable,dimension(:) :: source
real,allocatable,dimension(:) :: P1,P2,P3,VP
real,allocatable,dimension(:) :: Seism
integer :: index, index_src

! model parameters
Nx=1000
Nz=1000
h=10

inix=3
iniz=3
endx=Nx-2
endz=Nz-2

! time parameters
Nt=2001
dt=1.0e-3

!source parameters
sx   = Nx/2
sz   = 20
fcut = 30

! receiver depth
rz=20

! Allocate arrays
allocate(source(Nt))
allocate(VP(Nx*Nz))
allocate(P1(Nx*Nz))
allocate(P2(Nx*Nz))
allocate(P3(Nx*Nz))
allocate(Seism(Nt*Nx))


!Initializate arrays and counters
count_snap=1
Nsnap = 10
source = 0.0
P1 =0.0
P2 =0.0
P3 =0.0

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

        if (j<=200) then
            VP(index) = 1500
        else
            VP(index) = 1800
        end if

    end do
!$OMP END DO
!$OMP END PARALLEL

! Create Source Ricker
fcut_aux       = fcut/(3.*sqrt(pi))        ! Ajust to cut of gaussian function
t0_src   = 4*sqrt(pi)/fcut                 ! Initial time source
Nt_src = nint(2*t0_src/dt) + 1           ! Number of elements of the source
do k=1,Nt_src                          !Nts=nint(tm/dt)+1
    t_src=(k-1)*dt-t0_src                    !Delay Time
    src_aux=pi*(pi*fcut_aux*t_src)*(pi*fcut_aux*t_src)
    source(k) = (2*src_aux-1)*exp(-src_aux)    
end do

 ! Register in disk
 open(23, file='snapshots.bin', status='replace',&
 &FORM='unformatted',ACCESS='direct', recl=(Nx*Nz*4))

 open(24, file='seismogram.bin', status='replace',&
 &FORM='unformatted',ACCESS='direct', recl=(Nx*Nt*4))

 ! Solve wave equation
    do k=1,Nt
        index_src = sz + Nz*(sx-1)
        P2(index_src) = P2(index_src) + source(k)
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(P1,P2,P3,VP,Seism,dt,h,Nx,Nz,inix,iniz,endx,endz,Nt,rz)
        !$OMP DO
            ! Vectorized spatial loop    
            do index=1,Nx*Nz 
                !4th order in space and 2nd order in time
                if (mod(index,Nz)==0) then
                    i=index/Nz
                else    
                    i = index/Nz + 1
                end if

                j = index  - (i-1)*Nz

                if ((i>=inix .and. j>=iniz) .and. (i<endx .and. j< endz)) then            
                    ! P3(j + Nz*(i-1))=2*P2(j + Nz*(i-1))-P1(j + Nz*(i-1)) + ((VP(j + Nz*(i-1))*VP(j + Nz*(i-1)))*(dt*dt)/(12*h*h)) &
                    !         &*(-(P2(j + Nz*(i-2-1)) + P2(j-2 + Nz*(i-1)) + P2(j+2 + Nz*(i-1)) + P2(j + Nz*(i+2-1))) + & 
                    !         &16*(P2(j + Nz*(i-1-1)) + P2(j-1 + Nz*(i-1)) + P2(j+1 + Nz*(i-1)) + P2(j + Nz*(i+1-1))) - &
                    !         &60*P2(j + Nz*(i-1))) 
                    
                    p_zz = (-1*P2(j-2 + Nz*(i-1)) + 16*P2(j-1 + Nz*(i-1)) &
                    - 30*P2(j + Nz*(i-1)) + 16*P2(j+1 + Nz*(i-1)) - 1*P2(j+2 + Nz*(i-1) ) )/(12*h*h)

                    p_xx = (-1*P2(j + Nz*(i-1-2)) + 16*P2(j + Nz*(i-1-1)) &
                    - 30*P2(j + Nz*(i-1)) + 16*P2(j + Nz*(i-1+1)) - 1*P2(j + Nz*(i-1+2) ) ) / (12*h*h)
    
                    P3(index) = 2*P2(j + Nz*(i-1)) - P1(j + Nz*(i-1)) &
                    + (dt*dt)*(VP(j + Nz*(i-1))*VP(j + Nz*(i-1)))*(p_xx + p_zz)
                end if
            end do            
        !$OMP END DO NOWAIT
        
        ! update fields                                
        !$OMP DO
            do index=1,Nx*Nz
                P1(index)=P2(index)
                P2(index)=P3(index)    
            end do
        !$OMP END DO NOWAIT
        
        !Storage Seismogram        
        !$OMP DO
            do rx=1,Nx
                Seism(k + (rx-1)*Nt) = P3(rz + Nz*(rx-1))
            end do
        !$OMP END DO NOWAIT

        !$OMP END PARALLEL

        !Register snapshots
        !$OMP MASTER
            if (mod((k-1),Nt/Nsnap) == 0) then
                print*, "Propagation time =", (k-1)*dt, "Registering snapshot", count_snap
                write(23,rec=count_snap) ((P3(j + Nz*(i-1)) + 1.0e-3*VP(j + Nz*(i-1)),j=1,Nz),i=1,Nx)
                count_snap=count_snap+1
            end if     
        !$OMP END MASTER

       
    end do

    !Register Seismogram
    write(24,rec=1) (Seism(k),k=1,Nt*Nx)

!close files
close(23)
close(24)

end program