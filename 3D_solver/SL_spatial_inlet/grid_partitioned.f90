! ! Host-only main-domain mesh generation, independent of buffer parameters.
! module grid_partitioned
!   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
!   implicit none
!   private
!   public :: partitioned_axis
! contains
!   subroutine partitioned_axis(n, length, uniform_length, n_uniform, symmetric, coord)
!     integer, intent(in) :: n, n_uniform
!     real(8), intent(in) :: length, uniform_length
!     logical, intent(in) :: symmetric
!     real(8), intent(out) :: coord(n)
!     real(8) :: h, width(n-1), side_length
!     integer :: m, n_side, first, j
!     if (n < 2 .or. n_uniform < 2 .or. n_uniform > n) error stop 'Invalid uniform point count'
!     if (.not.ieee_is_finite(length) .or. .not.ieee_is_finite(uniform_length)) error stop 'Nonfinite mesh length'
!     if (uniform_length <= 0.d0 .or. uniform_length > length) error stop 'Invalid uniform length'
!     h = uniform_length / dble(n_uniform-1)
!     if (.not.symmetric) then
!       coord(1) = 0.d0
!       do j = 2, n_uniform
!         coord(j) = dble(j-1)*h
!       enddo
!       m = n-n_uniform
!       call stretched_widths(m, length-uniform_length, h, width(1:m))
!       do j = 1, m
!         coord(n_uniform+j) = coord(n_uniform+j-1)+width(j)
!       enddo
!     else
!       if (mod(n-n_uniform,2) /= 0) error stop 'Symmetric mesh requires even n-n_uniform'
!       n_side = (n-n_uniform)/2
!       first = n_side+1
!       side_length = 0.5d0*(length-uniform_length)
!       do j = 0, n_uniform-1
!         coord(first+j) = side_length+dble(j)*h
!       enddo
!       call stretched_widths(n_side, side_length, h, width(1:n_side))
!       do j = 1, n_side
!         coord(first-j) = coord(first-j+1)-width(j)
!         coord(first+n_uniform-1+j) = coord(first+n_uniform-2+j)+width(j)
!       enddo
!     endif
!     if (abs(coord(1)) > 1.d-11*length .or. abs(coord(n)-length) > 1.d-11*length) &
!       error stop 'Main mesh length mismatch'
!   end subroutine

!   subroutine stretched_widths(n, length, h, width)
!     integer, intent(in) :: n
!     real(8), intent(in) :: length, h
!     real(8), intent(out) :: width(n)
!     real(8) :: t, excess, weight(n), amplitude, tol
!     integer :: j
!     tol = 1.d-12*max(h, length)
!     excess = length-dble(n)*h
!     if (n == 0) then
!       if (abs(length) > tol) error stop 'Nonzero length without stretched intervals'
!       return
!     endif
!     if (excess < -tol) error stop 'Stretch section too short: increase uniform point count or reduce uniform length'
!     if (n == 1) then
!       if (abs(excess) > tol) error stop 'Need at least two stretched intervals for smooth matching'
!       width = h
!       return
!     endif
!     ! Width equals h at the junction; quintic ramp has zero slope at both ends.
!     do j = 1, n
!       t = dble(j-1)/dble(n-1)
!       weight(j) = t**3*(10.d0-15.d0*t+6.d0*t*t)
!     enddo
!     amplitude = max(0.d0, excess)/sum(weight)
!     width = h+amplitude*weight
!   end subroutine
! end module grid_partitioned
