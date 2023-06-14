! INPUT: REAL (single precision) datacube (n,n,n) where n is the number of grid per side
! n12 = n/2-1 for n is even ; IF n is odd, n12 = (n-1)/2 
! OUTPUT: D2power = dimensionless power spectrum k^3 P(k)/(2 pi^2) Vs unitless k = (1:n12). 
! physical k = (unitless k) * 2 pi/(boxsize) 
! In ranger, first run in command: module swap pgi intel; module load fftw3 -- to load mvapich and ifort and fftw3 
! compile like this: ifort -O3 -xW -openmp -I$TACC_FFTW3_INC nrtype.f90 powers.f90 mainbody.f90  -L$TACC_FFTW3_LIB -lfftw3_threads -lfftw3 -lm -o powers.x
! set OMP command in qsub command in Ranger


program mainbody
  use nrtype
  implicit none
  
  include '/home/andrewcaruso/fftw-3.3.9/include/fftw3.f'
  
  ! give n 
  ! if n is odd, change n12 here, and n121 in the subroutine compute_fftw
  integer, parameter :: n = 512.0, n12 = n/2-1 ! NORMALLY 512.0
  integer, parameter :: ind = 1 !1 = T21, !9 = gamma !16 = ion
  real(sp), parameter :: boxsize = 512.0
  real(sp), parameter :: kstep = 2.* PI / boxsize
  
  integer :: ierr, num_procs, OMP_GET_NUM_PROCS, kint
  character(len=200) :: dir, datafile, datafile1, datafile2, outfile
  real(sp),dimension(n,n,n) :: cube, cube1, cube2
  real(sp),dimension(1:n12) :: D2power 
  
  call dfftw_init_threads(ierr)
  if (ierr .eq. 0) stop 'Fatal error: multithread FFTW initiated incorrectly'
  num_procs = OMP_GET_NUM_PROCS()
  print *,'Number OMP PROCS = ', num_procs
  call dfftw_plan_with_nthreads(num_procs)

  ! give datafile here
  !dir = '/expanse/lustre/scratch/ngang002/temp_project/LAE_Bubble_Project_Scratch/Temperature21_Data/output_N3/'
  outfile = '/expanse/lustre/scratch/andrewcaruso/temp_project/final_powers/21cm_spectra_500subs/extreme/z=07.3402/cut_500/auto_21cm_=512:1024,0:512,512:1024'
  datafile1 = '/expanse/lustre/scratch/andrewcaruso/temp_project/final_powers/21cm_fields_500subs/extreme/z=07.3402/cut_500/t21_field.512:1024,0:512,512:1024'
  datafile2 = '/expanse/lustre/scratch/andrewcaruso/temp_project/final_powers/21cm_fields_500subs/extreme/z=07.3402/cut_500/t21_field.512:1024,0:512,512:1024'


  write(*,*) 'test'
  write(*,*) outfile

  ! auto power 
  !call read_density_cube2(n,datafile,cube)
  !call pipeline_w_spherical_averaging(n,n12,cube,D2power) 
  
  ! cross power of data 1 and data 2
  !call read_density_cube2(n,datafile1,cube1)
  !call read_density_cube2(n,datafile2,cube2)
  call read_gas_cube(n,datafile1,ind,cube1)
  write(*,*) 'read gas cube 1'
  call read_gas_cube(n,datafile2,ind,cube2)
  write(*,*) 'read gas cube 2'
  call pipeline_xpower_w_spherical_averaging(n,n12,cube1,cube2, D2power)
  
  ! then do whatever on D2power yourself
  open (20,file = outfile)
    write (20,'(a)')'k[unitless]    k[h/Mpc]    D^2(k)'
    do kint = 1,n12
      write (20,'(I6,20(E15.5))') kint, kstep*real(kint), D2power(kint)
    end do
  close (20)
  
  stop
end program mainbody
