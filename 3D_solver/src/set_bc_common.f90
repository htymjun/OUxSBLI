module set_bc_common
  implicit none
  interface set_bc_cyclic
    module procedure set_bc_cyclic2_init, set_bc_cyclic2, set_bc_cyclic4, &
                     set_bc_cyclic4_init, set_bc_cyclic6, set_bc_cyclic6_init
  end interface set_bc_cyclic
contains
  subroutine set_bc_cyclic2_init(id_accuracy, nx, ny, nz, Q)
    integer(kind=2), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout)             :: Q(5,nx,ny,nz)
    integer i, j, k
    do k = 2, nz-1
      do j = 2, ny-1
        Q(:,1,j,k)  = Q(:,nx-1,j,k)
        Q(:,nx,j,k) = Q(:,2,j,k)
    enddo;enddo
    do k = 2, nz-1
      do i = 2, nx-1
        Q(:,i,1,k)  = Q(:,i,ny-1,k)
        Q(:,i,ny,k) = Q(:,i,2,k)
    enddo;enddo
    do k = 2, nz-1
      Q(:,1,1,k)   = Q(:,nx-1,ny-1,k)
      Q(:,nx,1,k)  = Q(:,2,ny-1,k)
      Q(:,1,ny,k)  = Q(:,nx-1,2,k)
      Q(:,nx,ny,k) = Q(:,2,2,k)
    enddo
    do j = 1, ny
      do i = 1, nx
        Q(:,i,j,1)  = Q(:,i,j,nz-1)
        Q(:,i,j,nz) = Q(:,i,j,2)
    enddo;enddo
  end subroutine set_bc_cyclic2_init 


  subroutine set_bc_cyclic4_init(id_accuracy, nx, ny, nz, Q)
    integer(kind=4), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout)             :: Q(5,nx,ny,nz)
    integer i, j, k
    do k = 3, nz-2
      do j = 3, ny-2
        Q(:,1:2,j,k) = Q(:,nx-3:nx-2,j,k)
        Q(:,nx-1:nx,j,k) = Q(:,3:4,j,k)
    enddo;enddo
    do k = 3, nz-2
      do i = 3, nx-2
        Q(:,i,1:2,k) = Q(:,i,ny-3:ny-2,k)
        Q(:,i,ny-1:ny,k) = Q(:,i,3:4,k)
    enddo;enddo
    do k = 3, nz-2
      Q(:,1:2,1:2,k) = Q(:,nx-3:nx-2,ny-3:ny-2,k)
      Q(:,nx-1:nx,1:2,k) = Q(:,3:4,ny-3:ny-2,k)
      Q(:,1:2,ny-1:ny,k) = Q(:,nx-3:nx-2,3:4,k)
      Q(:,nx-1:nx,ny-1:ny,k) = Q(:,3:4,3:4,k)
    enddo
    do j = 1, ny
      do i = 1, nx
        Q(:,i,j,1:2) = Q(:,i,j,nz-3:nz-2)
        Q(:,i,j,nz-1:nz) = Q(:,i,j,3:4)
    enddo;enddo
  end subroutine set_bc_cyclic4_init


  subroutine set_bc_cyclic6_init(id_accuracy, nx, ny, nz, Q)
    integer(kind=8), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout)             :: Q(5,nx,ny,nz)
    integer i, j, k
    do k = 4, nz-3
      do j = 4, ny-3
        Q(:,1:3,j,k) = Q(:,nx-5:nx-3,j,k)
        Q(:,nx-2:nx,j,k) = Q(:,4:6,j,k)
    enddo;enddo
    do k = 4, nz-3
      do i = 4, nx-3
        Q(:,i,1:3,k) = Q(:,i,ny-5:ny-3,k)
        Q(:,i,ny-2:ny,k) = Q(:,i,4:6,k)
    enddo;enddo
    do k = 4, nz-3
      Q(:,1:3,1:3,k) = Q(:,nx-5:nx-3,ny-5:ny-3,k)
      Q(:,nx-2:nx,1:3,k) = Q(:,4:6,ny-5:ny-3,k)
      Q(:,1:3,ny-2:ny,k) = Q(:,nx-5:nx-3,4:6,k)
      Q(:,nx-2:nx,ny-2:ny,k) = Q(:,4:6,4:6,k)
    enddo
    do j = 1, ny
      do i = 1, nx
        Q(:,i,j,1:3) = Q(:,i,j,nz-5:nz-3)
        Q(:,i,j,nz-2:nz) = Q(:,i,j,4:6)
    enddo;enddo
  end subroutine set_bc_cyclic6_init


  subroutine set_bc_cyclic2(id_accuracy, nx, ny, nz, Q)
    integer(kind=2), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout), device     :: Q(5,nx,ny,nz)
    integer i, j, k, l
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        do l = 1, 5
          Q(l,1,j,k)  = Q(l,nx-1,j,k)
          Q(l,nx,j,k) = Q(l,2,j,k)
    enddo;enddo;enddo
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do i = 2, nx-1
        do l = 1, 5
          Q(l,i,1,k)  = Q(l,i,ny-1,k)
          Q(l,i,ny,k) = Q(l,i,2,k)
    enddo;enddo;enddo
    !$cuf kernel do(1)<<<*,*>>>
    do k = 2, nz-1
      do l = 1, 5
        Q(l,1,1,k)   = Q(l,nx-1,ny-1,k)
        Q(l,nx,1,k)  = Q(l,2,ny-1,k)
        Q(l,1,ny,k)  = Q(l,nx-1,2,k)
        Q(l,nx,ny,k) = Q(l,2,2,k)
    enddo;enddo
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        do l = 1, 5
          Q(l,i,j,1)  = Q(l,i,j,nz-1)
          Q(l,i,j,nz) = Q(l,i,j,2)
    enddo;enddo;enddo
  end subroutine set_bc_cyclic2


  subroutine set_bc_cyclic4(id_accuracy, nx, ny, nz, Q)
    integer(kind=4), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout), device     :: Q(5,nx,ny,nz)
    integer i, j, k, l
    !$cuf kernel do(2) <<<*,*>>>
    do k = 3, nz-2
      do j = 3, ny-2
        do l = 1, 5
          Q(l,1,j,k) = Q(l,nx-3,j,k)
          Q(l,2,j,k) = Q(l,nx-2,j,k)
          Q(l,nx-1,j,k) = Q(l,3,j,k)
          Q(l,nx,j,k) = Q(l,4,j,k)
    enddo;enddo;enddo
    !$cuf kernel do(2) <<<*,*>>>
    do k = 3, nz-2
      do i = 3, nx-2
        do l = 1, 5
          Q(l,i,1,k) = Q(l,i,ny-3,k)
          Q(l,i,2,k) = Q(l,i,ny-2,k)
          Q(l,i,ny-1,k) = Q(l,i,3,k)
          Q(l,i,ny,k) = Q(l,i,4,k)
    enddo;enddo;enddo
    !$cuf kernel do(1) <<<*,*>>>
    do k = 3, nz-2
      do l = 1, 5
        Q(l,1,1,k) = Q(l,nx-3,ny-3,k)
        Q(l,1,2,k) = Q(l,nx-3,ny-2,k)
        Q(l,2,1,k) = Q(l,nx-2,ny-3,k)
        Q(l,2,2,k) = Q(l,nx-2,ny-2,k)
        Q(l,nx-1,1,k) = Q(l,3,ny-3,k)
        Q(l,nx-1,2,k) = Q(l,3,ny-2,k)
        Q(l,nx,1,k)   = Q(l,4,ny-3,k)
        Q(l,nx,2,k)   = Q(l,4,ny-2,k)
        Q(l,1,ny-1,k) = Q(l,nx-3,3,k)
        Q(l,1,ny,k)   = Q(l,nx-3,4,k)
        Q(l,2,ny-1,k) = Q(l,nx-2,3,k)
        Q(l,2,ny,k)   = Q(l,nx-2,4,k)
        Q(l,nx-1,ny-1,k) = Q(l,3,3,k)
        Q(l,nx-1,ny,k)   = Q(l,3,4,k)
        Q(l,nx,ny-1,k)   = Q(l,4,3,k)
        Q(l,nx,ny,k)     = Q(l,4,4,k)
    enddo;enddo
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        do l = 1, 5
          Q(l,i,j,1) = Q(l,i,j,nz-3)
          Q(l,i,j,2) = Q(l,i,j,nz-2)
          Q(l,i,j,nz-1) = Q(l,i,j,3)
          Q(l,i,j,nz) = Q(l,i,j,4)
    enddo;enddo;enddo
  end subroutine set_bc_cyclic4


  subroutine set_bc_cyclic6(id_accuracy, nx, ny, nz, Q)
    integer(kind=8), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout), device     :: Q(5,nx,ny,nz)
    integer i, j, k, l
    !$cuf kernel do(2)<<<*,*>>>
    do k = 4, nz-3
      do j = 4, ny-3
        do l = 1, 5
          Q(l,1,j,k) = Q(l,nx-5,j,k)
          Q(l,2,j,k) = Q(l,nx-4,j,k)
          Q(l,3,j,k) = Q(l,nx-3,j,k)
          Q(l,nx-2,j,k) = Q(l,4,j,k)
          Q(l,nx-1,j,k) = Q(l,5,j,k)
          Q(l,nx,j,k)   = Q(l,6,j,k)
    enddo;enddo;enddo
    !$cuf kernel do(2)<<<*,*>>>
    do k = 4, nz-3
      do i = 4, nx-3
        do l = 1, 5
          Q(l,i,1,k) = Q(l,i,ny-5,k)
          Q(l,i,2,k) = Q(l,i,ny-4,k)
          Q(l,i,3,k) = Q(l,i,ny-3,k)
          Q(l,i,ny-2,k) = Q(l,i,4,k)
          Q(l,i,ny-1,k) = Q(l,i,5,k)
          Q(l,i,ny,k)   = Q(l,i,6,k)
    enddo;enddo;enddo
    !$cuf kernel do(1)<<<*,*>>>
    do k = 4, nz-3
      do l = 1, 5
        Q(l,1,1,k) = Q(l,nx-5,ny-5,k)
        Q(l,1,2,k) = Q(l,nx-5,ny-4,k)
        Q(l,1,3,k) = Q(l,nx-5,ny-3,k)
        Q(l,2,1,k) = Q(l,nx-4,ny-5,k)
        Q(l,2,2,k) = Q(l,nx-4,ny-4,k)
        Q(l,2,3,k) = Q(l,nx-4,ny-3,k)
        Q(l,3,1,k) = Q(l,nx-3,ny-5,k)
        Q(l,3,2,k) = Q(l,nx-3,ny-4,k)
        Q(l,3,3,k) = Q(l,nx-3,ny-3,k)
        Q(l,nx-2,1,k) = Q(l,4,ny-5,k)
        Q(l,nx-2,2,k) = Q(l,4,ny-4,k)
        Q(l,nx-2,3,k) = Q(l,4,ny-3,k)
        Q(l,nx-1,1,k) = Q(l,5,ny-5,k)
        Q(l,nx-1,2,k) = Q(l,5,ny-4,k)
        Q(l,nx-1,3,k) = Q(l,5,ny-3,k)
        Q(l,nx,1,k)   = Q(l,6,ny-5,k)
        Q(l,nx,2,k)   = Q(l,6,ny-4,k)
        Q(l,nx,3,k)   = Q(l,6,ny-3,k)
        Q(l,1,ny-2,k) = Q(l,nx-5,4,k)
        Q(l,1,ny-1,k) = Q(l,nx-5,5,k)
        Q(l,1,ny,k)   = Q(l,nx-5,6,k)
        Q(l,2,ny-2,k) = Q(l,nx-4,4,k)
        Q(l,2,ny-1,k) = Q(l,nx-4,5,k)
        Q(l,2,ny,k)   = Q(l,nx-4,6,k)
        Q(l,3,ny-2,k) = Q(l,nx-3,4,k)
        Q(l,3,ny-1,k) = Q(l,nx-3,5,k)
        Q(l,3,ny,k)   = Q(l,nx-3,6,k)
        Q(l,nx-2,ny-2,k) = Q(l,4,4,k)
        Q(l,nx-2,ny-1,k) = Q(l,4,5,k)
        Q(l,nx-2,ny,k)   = Q(l,4,6,k)
        Q(l,nx-1,ny-2,k) = Q(l,5,4,k)
        Q(l,nx-1,ny-1,k) = Q(l,5,5,k)
        Q(l,nx-1,ny,k)   = Q(l,5,6,k)
        Q(l,nx,ny-2,k)   = Q(l,6,4,k)
        Q(l,nx,ny-1,k)   = Q(l,6,5,k)
        Q(l,nx,ny,k)     = Q(l,6,6,k)
    enddo;enddo
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        do l = 1, 5
          Q(l,i,j,1) = Q(l,i,j,nz-5)
          Q(l,i,j,2) = Q(l,i,j,nz-4)
          Q(l,i,j,3) = Q(l,i,j,nz-3)
          Q(l,i,j,nz-2) = Q(l,i,j,4)
          Q(l,i,j,nz-1) = Q(l,i,j,5)
          Q(l,i,j,nz)   = Q(l,i,j,6)
    enddo;enddo;enddo
  end subroutine set_bc_cyclic6


  subroutine set_bc_cyclic_z(nx, ny, nz, QJ)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: QJ(5,nx,ny,nz)
    integer i, j, l
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        do l = 1, 5
          QJ(l,i,j,1) = QJ(l,i,j,nz-5)
          QJ(l,i,j,2) = QJ(l,i,j,nz-4)
          QJ(l,i,j,3) = QJ(l,i,j,nz-3)
          QJ(l,i,j,nz-2) = QJ(l,i,j,4)
          QJ(l,i,j,nz-1) = QJ(l,i,j,5)
          QJ(l,i,j,nz)   = QJ(l,i,j,6)
    enddo;enddo;enddo
  end subroutine set_bc_cyclic_z


  subroutine set_bc_mut_common(nx, ny, nz, mut, qc2)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: mut(nx,ny,nz), qc2(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(2) <<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        mut(1,j,k) = mut(nx-1,j,k)
        mut(nx,j,k) = mut(2,j,k)
        qc2(1,j,k) = qc2(nx-1,j,k)
        qc2(nx,j,k) = qc2(2,j,k)
    enddo;enddo
    !$cuf kernel do(2) <<<*,*>>>
    do k = 2, nz-1
      do i = 2, nx-1
        mut(i,1,k) = mut(i,ny-1,k)
        mut(i,ny,k) = mut(i,2,k)
        qc2(i,1,k) = qc2(i,ny-1,k)
        qc2(i,ny,k) = qc2(i,2,k)
    enddo;enddo
    !$cuf kernel do(1) <<<*,*>>>
    do k = 2, nz-1
      mut(1,1,k) = mut(nx-1,ny-1,k)
      mut(nx,1,k) = mut(2,ny-1,k)
      mut(1,ny,k) = mut(nx-1,2,k)
      mut(nx,ny,k) = mut(2,2,k)
      qc2(1,1,k) = qc2(nx-1,ny-1,k)
      qc2(nx,1,k) = qc2(2,ny-1,k)
      qc2(1,ny,k) = qc2(nx-1,2,k)
      qc2(nx,ny,k) = qc2(2,2,k)
    enddo
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        mut(i,j,1) = mut(i,j,nz-1)
        mut(i,j,nz) = mut(i,j,2)
        qc2(i,j,1) = qc2(i,j,nz-1)
        qc2(i,j,nz) = qc2(i,j,2)
    enddo;enddo
  end subroutine set_bc_mut_common
end module set_bc_common

