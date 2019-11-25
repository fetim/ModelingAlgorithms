Program VEC_ADD_DO_OPENMP
!$ use omp_lib
implicit none

integer, parameter :: n=10
integer :: i, a(n),b(n),t(n),tid,stack
!$OMP PARALLEL SHARED(a,t) PRIVATE(i,tid) REDUCTION(+:stack)
    tid = 0
    a=0
    b=0
    stack =0 
    !$ tid = OMP_GET_THREAD_NUM()
    !$OMP DO    
        do i=1,n
            a(i) = i
            b(i) = stack
            stack = stack + a(i)
            t(i) = tid ! Record whic thread did which iteration
        end do
    !$OMP END DO
!$OMP END PARALLEL

    do i =1,n
        write(*,*) i, a(i), b(i), t(i)
        
    end do
    
End program