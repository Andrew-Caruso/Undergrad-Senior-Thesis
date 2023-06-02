! Yi Mao 2008-2010

subroutine read_density_cube2(n,datafile,cube) !rho, not delta_rho
! read number density 
! file format: header (no extras), single precision
  use nrtype
  implicit none
  integer, intent(in) :: n
  character(len=200),intent(in) :: datafile
  real(sp) :: time
  real(sp),dimension(n,n,n),intent(out) :: cube 
  integer :: ndread

  open (10,form='binary',file=datafile)
    !read (10) ndread, ndread, ndread, cube
    read (10) cube
  close (10)
  !if (ndread .ne. n) stop 'Reading Fator error: dimension not matching'
  return
end subroutine read_density_cube2

subroutine read_gas_cube(n,datafile,ind,cube) !rho, not delta_rho
! read number density 
! file format: header (no extras), single precision
  use nrtype
  implicit none
  integer, intent(in) :: n, ind
  character(len=200),intent(in) :: datafile
  real(sp) :: time
  real(sp),dimension(1,n,n,n) :: gas_cube
  real(sp),dimension(n,n,n),intent(out) :: cube
  integer :: ndread

  open (10,form='binary',file=datafile)
    !read (10) ndread, ndread, ndread, cube
    !read (10) time
    read (10) gas_cube
  close (10)
  !write(*,*) time
  write(*,*) 'test 1'
  cube = gas_cube(ind, :, :, :)
  write(*,*) cube(1, 1, 1)
  !if (ndread .ne. n) stop 'Reading Fator error: dimension not matching'
  return
end subroutine read_gas_cube

!--------- pipelines --------------------
! n12 = n/2-1 for n is even ; n12 = (n-1)/2  for n is odd
subroutine pipeline_w_spherical_averaging(n,n12,cube,D2power)
! compute auto power spectrm 
  use nrtype
  implicit none
  integer, intent(in) :: n, n12
  real(sp),dimension(n,n,n),intent(in):: cube 
  real(sp),dimension(1:n12),intent(out) :: D2power 
  double complex, dimension(:,:,:),allocatable :: cubeF
  real(sp),dimension(:,:,:),allocatable:: PowerW
  real(sp),dimension(:),allocatable ::  dev
  integer, dimension(:),allocatable :: Nmodes
  
  allocate(cubeF(0:n12,-n12:n12,-n12:n12),PowerW(0:n12,-n12:n12,-n12:n12)  & 
    ,dev(1:n12),Nmodes(1:n12))
  call compute_fftw(n,n12,cube,cubeF)
  call compute_single_power_spectrum(n12,cubeF,cubeF,PowerW)
  call simple_spherical_averaging_over_modes(n,n12,PowerW,D2power,dev,Nmodes)
  deallocate(cubeF,PowerW,dev,Nmodes)
  return
end subroutine pipeline_w_spherical_averaging


subroutine pipeline_xpower_w_spherical_averaging(n,n12,cube1,cube2, D2power)
! compute cross power spectrum P12 = 1/2*( cc(d1) * d2 + d1 * cc(d2)) where cc = complex conjugate
  use nrtype
  implicit none
  integer, intent(in) :: n, n12
  real(sp),dimension(n,n,n),intent(in):: cube1, cube2
  real(sp),dimension(1:n12),intent(out) :: D2power 
  double complex, dimension(:,:,:),allocatable :: cube1F, cube2F
  real(sp),dimension(:,:,:),allocatable:: PowerW
  real(sp),dimension(:),allocatable ::  dev
  integer, dimension(:),allocatable :: Nmodes
  
  allocate(cube1F(0:n12,-n12:n12,-n12:n12),cube2F(0:n12,-n12:n12,-n12:n12) & 
  ,PowerW(0:n12,-n12:n12,-n12:n12),dev(1:n12),Nmodes(1:n12))
  call compute_fftw(n,n12,cube1,cube1F)
  call compute_fftw(n,n12,cube2,cube2F)
  call compute_single_power_spectrum(n12,cube1F,cube2F,PowerW)
  call simple_spherical_averaging_over_modes(n,n12,PowerW,D2power,dev,Nmodes)
  deallocate(cube1F,cube2F,PowerW,dev,Nmodes)
  return
end subroutine pipeline_xpower_w_spherical_averaging

!------------- compute the FFT of a data cube (unnormalized FFTW) -----------------
subroutine compute_fftw(n,n12,in,out)
! n12 = n/2-1 for n is even ; n12 = (n-1)/2  for n is odd
! note: if n = even, modes at the outest boundaries k=nd121 are abandoned.
  use nrtype
  implicit none

  include '/home/ngang002/fftw-3.3.9/include/fftw3.f'
  
  integer, intent(in):: n, n12
  real(sp),dimension(1:n,1:n,1:n),intent(in):: in
  double complex,dimension(:,:,:),allocatable:: inter
  double complex,dimension(0:n12,-n12:n12,-n12:n12),intent(out) :: out
  ! for n=even, some modes at n121 are not in the output cube
  integer*8:: plan
  integer :: n121, nm1, j1,j2,j3,i1,i2,i3
  real(dp),dimension(:,:,:),allocatable :: indb !(1:n,1:n,1:n)
  
  nm1 = n - 1
  n121 = n12+1 !!!nd is even 
  !!!n121= n12   !nd is odd
  
  allocate(indb(1:n,1:n,1:n),inter(0:n121,0:nm1,0:nm1))
  indb = dble(in)
  call dfftw_plan_dft_r2c_3d(plan,n,n,n,indb,inter,FFTW_ESTIMATE)
  call dfftw_execute_dft_r2c(plan, indb, inter)
  call dfftw_destroy_plan(plan)
  deallocate(indb)

  !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(n12,n,out,inter) PRIVATE(j3,j2,j1,i3,i2,i1)
  do j3 = -n12, n12
    do j2 = -n12, n12
      do j1 = 0,n12
        i1 = j1
        if (j2 .lt. 0 ) then 
          i2 = j2+n
        else 
          i2 = j2
        endif

        if (j3 .lt. 0 ) then
          i3 = j3+n
        else
          i3 = j3
        endif
        out(j1,j2,j3) = inter(i1,i2,i3)
      end do
    end do  
  end do
  !$OMP END PARALLEL DO
  deallocate(inter)
  return
end subroutine compute_fftw


!------------- compute (unnormalized) power spectra ------------------
subroutine compute_single_power_spectrum(n12,fF,gF,PkW)
  use nrtype
  implicit none
  integer, intent(in) :: n12
  double complex,dimension(0:n12,-n12:n12,-n12:n12),intent(in) :: fF,gF
  real(sp),dimension(0:n12,-n12:n12,-n12:n12),intent(out) :: PkW
  integer :: k1,k2,k3
  
  !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(n12,fF,gF,PkW) PRIVATE(k3,k2,k1)
  do k3 = -n12, n12
    do k2 = -n12, n12
      do k1 = 0, n12  
        PkW(k1,k2,k3) = real(CONJG(fF(k1,k2,k3))*gF(k1,k2,k3)) 
      end do
    end do
  end do
  !$OMP END PARALLEL DO
  return
end subroutine compute_single_power_spectrum

!----------- Power spectrum averaging ------------------
subroutine simple_spherical_averaging_over_modes(n,n12,cube,D2power,dev,Nmodes)
! cube is power spectrum in a half-cube in FFTW unit convention
! D2power is Delta^2 (unitless mean power in an annulus); dev=(stand.dev.)/(mean); count=#modes in a shell
  use nrtype
  implicit none
  integer, intent(in) :: n,n12
  real(sp),dimension(0:n12,-n12:n12,-n12:n12),intent(in) :: cube
  real(sp),dimension(1:n12),intent(out) :: D2power, dev
  integer, dimension(1:n12),intent(out) :: Nmodes
  real(dp),dimension(1:n12) :: power
  real(dp),dimension(1:n12) :: powersum, sqsum,sq
  integer :: k1,k2,k3, k

  powersum = 0.d0
  sqsum = 0.d0
  Nmodes = 0
  do k3 = -n12, n12
    do k2 = -n12, n12
      do k1 = 0, n12  
        k = NINT(sqrt(dble(k1)**2.+dble(k2)**2.+dble(k3)**2.))  !find nearest integer
        if (k .ge. 1 .and. k .le. n12) then
          if (k1 .ne. 0) then
            Nmodes(k) = Nmodes(k) + 2
            powersum(k) = powersum(k) + dble(cube(k1,k2,k3)) * 2.
            sqsum(k) = sqsum(k) + (dble(cube(k1,k2,k3)))**2. * 2.
          else 
            Nmodes(k) = Nmodes(k) + 1
            powersum(k) = powersum(k) + dble(cube(k1,k2,k3))
            sqsum(k) = sqsum(k) + (dble(cube(k1,k2,k3)))**2.
          endif 
        endif
      end do
    end do
  end do
  power = powersum/dble(Nmodes) 
  sq = sqsum/dble(Nmodes)
  dev = real(sqrt(sq/power**2. - 1.))
  !$OMP PARALLEL DO SCHEDULE(STATIC) SHARED(n12,n,power,D2power) PRIVATE(k)
  do k = 1, n12 
    D2power(k) = real(power(k) * 4.*pi_d/dble(n)**6. * dble(k)**3.)
  end do
  !$OMP END PARALLEL DO
  return
end subroutine simple_spherical_averaging_over_modes
