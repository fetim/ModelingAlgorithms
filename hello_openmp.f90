Program helloopenmp
implicit none

    integer nthreads, tid
    integer OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

    ! Fork a team of threads giving them their own copies of variables
    !$OMP PARALLEL PRIVATE(NTHREADS,tid)

    ! Obtain and print thread id
    !$ tid = OMP_GET_THREAD_NUM()
    print*, "Hello world from thread =", tid

    ! Only master thread does this
    if (tid == 0) then
       !$ NTHREADS = OMP_GET_NUM_THREADS()
        print*, "Number of threads=", NTHREADS
    end if
    !$OMP END PARALLEL
end program 