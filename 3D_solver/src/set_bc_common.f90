module set_bc_common
  use mod_globals, only : nx, ny, nz
  implicit none
  interface set_bc_cyclic
    module procedure set_bc_cyclic2_init, set_bc_cyclic2, set_bc_cyclic4, &
                     set_bc_cyclic4_init, set_bc_cyclic6, set_bc_cyclic6_init
  end interface set_bc_cyclic
contains
  !> Cyclic boundary condition initialization for second-order accuracy
  !> Executed on CPU before main time-stepping loop
  subroutine set_bc_cyclic2_init(id_accuracy, nx, ny, nz, Q)
    integer(kind=2), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout)             :: Q(nx,ny,nz,5)
    integer i, j, k
    do k = 2, nz-1
      do j = 2, ny-1
        Q(1,j,k,:)  = Q(nx-1,j,k,:)
        Q(nx,j,k,:) = Q(2,j,k,:)
    enddo;enddo
    do k = 2, nz-1
      do i = 2, nx-1
        Q(i,1,k,:)  = Q(i,ny-1,k,:)
        Q(i,ny,k,:) = Q(i,2,k,:)
    enddo;enddo
    do k = 2, nz-1
      Q(1,1,k,:)   = Q(nx-1,ny-1,k,:)
      Q(nx,1,k,:)  = Q(2,ny-1,k,:)
      Q(1,ny,k,:)  = Q(nx-1,2,k,:)
      Q(nx,ny,k,:) = Q(2,2,k,:)
    enddo
    do j = 1, ny
      do i = 1, nx
        Q(i,j,1,:)  = Q(i,j,nz-1,:)
        Q(i,j,nz,:) = Q(i,j,2,:)
    enddo;enddo
  end subroutine set_bc_cyclic2_init

  !> Cyclic boundary condition initialization for 4th-order accuracy
  !> Executed on CPU before main time-stepping loop
  subroutine set_bc_cyclic4_init(id_accuracy, nx, ny, nz, Q)
    integer(kind=4), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout)             :: Q(nx,ny,nz,5)
    integer i, j, k
    do k = 3, nz-2
      do j = 3, ny-2
        Q(1:2,j,k,:) = Q(nx-3:nx-2,j,k,:)
        Q(nx-1:nx,j,k,:) = Q(3:4,j,k,:)
    enddo;enddo
    do k = 3, nz-2
      do i = 3, nx-2
        Q(i,1:2,k,:) = Q(i,ny-3:ny-2,k,:)
        Q(i,ny-1:ny,k,:) = Q(i,3:4,k,:)
    enddo;enddo
    do k = 3, nz-2
      Q(1:2,1:2,k,:) = Q(nx-3:nx-2,ny-3:ny-2,k,:)
      Q(nx-1:nx,1:2,k,:) = Q(3:4,ny-3:ny-2,k,:)
      Q(1:2,ny-1:ny,k,:) = Q(nx-3:nx-2,3:4,k,:)
      Q(nx-1:nx,ny-1:ny,k,:) = Q(3:4,3:4,k,:)
    enddo
    do j = 1, ny
      do i = 1, nx
        Q(i,j,1:2,:) = Q(i,j,nz-3:nz-2,:)
        Q(i,j,nz-1:nz,:) = Q(i,j,3:4,:)
    enddo;enddo
  end subroutine set_bc_cyclic4_init

  !> Cyclic boundary condition initialization for 6th-order accuracy
  !> Executed on CPU before main time-stepping loop
  subroutine set_bc_cyclic6_init(id_accuracy, nx, ny, nz, Q)
    integer(kind=8), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout)             :: Q(nx,ny,nz,5)
    integer i, j, k
    do k = 4, nz-3
      do j = 4, ny-3
        Q(1:3,j,k,:) = Q(nx-5:nx-3,j,k,:)
        Q(nx-2:nx,j,k,:) = Q(4:6,j,k,:)
    enddo;enddo
    do k = 4, nz-3
      do i = 4, nx-3
        Q(i,1:3,k,:) = Q(i,ny-5:ny-3,k,:)
        Q(i,ny-2:ny,k,:) = Q(i,4:6,k,:)
    enddo;enddo
    do k = 4, nz-3
      Q(1:3,1:3,k,:) = Q(nx-5:nx-3,ny-5:ny-3,k,:)
      Q(nx-2:nx,1:3,k,:) = Q(4:6,ny-5:ny-3,k,:)
      Q(1:3,ny-2:ny,k,:) = Q(nx-5:nx-3,4:6,k,:)
      Q(nx-2:nx,ny-2:ny,k,:) = Q(4:6,4:6,k,:)
    enddo
    do j = 1, ny
      do i = 1, nx
        Q(i,j,1:3,:) = Q(i,j,nz-5:nz-3,:)
        Q(i,j,nz-2:nz,:) = Q(i,j,4:6,:)
    enddo;enddo
  end subroutine set_bc_cyclic6_init

  !> Cyclic boundary condition for second-order accuracy
  !> Executed on GPU during main time-stepping loop
  subroutine set_bc_cyclic2(id_accuracy, nx, ny, nz, Q_1, Q_2, Q_3, Q_4, Q_5)
    integer(kind=2), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout), device     :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do j = 2, ny-1
        Q_1(1,j,k)  = Q_1(nx-1,j,k);  Q_1(nx,j,k) = Q_1(2,j,k)
        Q_2(1,j,k)  = Q_2(nx-1,j,k);  Q_2(nx,j,k) = Q_2(2,j,k)
        Q_3(1,j,k)  = Q_3(nx-1,j,k);  Q_3(nx,j,k) = Q_3(2,j,k)
        Q_4(1,j,k)  = Q_4(nx-1,j,k);  Q_4(nx,j,k) = Q_4(2,j,k)
        Q_5(1,j,k)  = Q_5(nx-1,j,k);  Q_5(nx,j,k) = Q_5(2,j,k)
    enddo;enddo
    !$cuf kernel do(2)<<<*,*>>>
    do k = 2, nz-1
      do i = 2, nx-1
        Q_1(i,1,k)  = Q_1(i,ny-1,k);  Q_1(i,ny,k) = Q_1(i,2,k)
        Q_2(i,1,k)  = Q_2(i,ny-1,k);  Q_2(i,ny,k) = Q_2(i,2,k)
        Q_3(i,1,k)  = Q_3(i,ny-1,k);  Q_3(i,ny,k) = Q_3(i,2,k)
        Q_4(i,1,k)  = Q_4(i,ny-1,k);  Q_4(i,ny,k) = Q_4(i,2,k)
        Q_5(i,1,k)  = Q_5(i,ny-1,k);  Q_5(i,ny,k) = Q_5(i,2,k)
    enddo;enddo
    !$cuf kernel do(1)<<<*,*>>>
    do k = 2, nz-1
      Q_1(1,1,k)   = Q_1(nx-1,ny-1,k);  Q_1(nx,1,k)  = Q_1(2,ny-1,k)
      Q_1(1,ny,k)  = Q_1(nx-1,2,k);     Q_1(nx,ny,k) = Q_1(2,2,k)
      Q_2(1,1,k)   = Q_2(nx-1,ny-1,k);  Q_2(nx,1,k)  = Q_2(2,ny-1,k)
      Q_2(1,ny,k)  = Q_2(nx-1,2,k);     Q_2(nx,ny,k) = Q_2(2,2,k)
      Q_3(1,1,k)   = Q_3(nx-1,ny-1,k);  Q_3(nx,1,k)  = Q_3(2,ny-1,k)
      Q_3(1,ny,k)  = Q_3(nx-1,2,k);     Q_3(nx,ny,k) = Q_3(2,2,k)
      Q_4(1,1,k)   = Q_4(nx-1,ny-1,k);  Q_4(nx,1,k)  = Q_4(2,ny-1,k)
      Q_4(1,ny,k)  = Q_4(nx-1,2,k);     Q_4(nx,ny,k) = Q_4(2,2,k)
      Q_5(1,1,k)   = Q_5(nx-1,ny-1,k);  Q_5(nx,1,k)  = Q_5(2,ny-1,k)
      Q_5(1,ny,k)  = Q_5(nx-1,2,k);     Q_5(nx,ny,k) = Q_5(2,2,k)
    enddo
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        Q_1(i,j,1)  = Q_1(i,j,nz-1);  Q_1(i,j,nz) = Q_1(i,j,2)
        Q_2(i,j,1)  = Q_2(i,j,nz-1);  Q_2(i,j,nz) = Q_2(i,j,2)
        Q_3(i,j,1)  = Q_3(i,j,nz-1);  Q_3(i,j,nz) = Q_3(i,j,2)
        Q_4(i,j,1)  = Q_4(i,j,nz-1);  Q_4(i,j,nz) = Q_4(i,j,2)
        Q_5(i,j,1)  = Q_5(i,j,nz-1);  Q_5(i,j,nz) = Q_5(i,j,2)
    enddo;enddo
  end subroutine set_bc_cyclic2

  !> Cyclic boundary condition for 4th-order accuracy
  !> Executed on GPU during main time-stepping loop
  subroutine set_bc_cyclic4(id_accuracy, nx, ny, nz, Q_1, Q_2, Q_3, Q_4, Q_5)
    integer(kind=4), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout), device     :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(2) <<<*,*>>>
    do k = 3, nz-2
      do j = 3, ny-2
        Q_1(1,j,k) = Q_1(nx-3,j,k); Q_1(2,j,k) = Q_1(nx-2,j,k); Q_1(nx-1,j,k) = Q_1(3,j,k); Q_1(nx,j,k) = Q_1(4,j,k)
        Q_2(1,j,k) = Q_2(nx-3,j,k); Q_2(2,j,k) = Q_2(nx-2,j,k); Q_2(nx-1,j,k) = Q_2(3,j,k); Q_2(nx,j,k) = Q_2(4,j,k)
        Q_3(1,j,k) = Q_3(nx-3,j,k); Q_3(2,j,k) = Q_3(nx-2,j,k); Q_3(nx-1,j,k) = Q_3(3,j,k); Q_3(nx,j,k) = Q_3(4,j,k)
        Q_4(1,j,k) = Q_4(nx-3,j,k); Q_4(2,j,k) = Q_4(nx-2,j,k); Q_4(nx-1,j,k) = Q_4(3,j,k); Q_4(nx,j,k) = Q_4(4,j,k)
        Q_5(1,j,k) = Q_5(nx-3,j,k); Q_5(2,j,k) = Q_5(nx-2,j,k); Q_5(nx-1,j,k) = Q_5(3,j,k); Q_5(nx,j,k) = Q_5(4,j,k)
    enddo;enddo
    !$cuf kernel do(2) <<<*,*>>>
    do k = 3, nz-2
      do i = 3, nx-2
        Q_1(i,1,k) = Q_1(i,ny-3,k); Q_1(i,2,k) = Q_1(i,ny-2,k); Q_1(i,ny-1,k) = Q_1(i,3,k); Q_1(i,ny,k) = Q_1(i,4,k)
        Q_2(i,1,k) = Q_2(i,ny-3,k); Q_2(i,2,k) = Q_2(i,ny-2,k); Q_2(i,ny-1,k) = Q_2(i,3,k); Q_2(i,ny,k) = Q_2(i,4,k)
        Q_3(i,1,k) = Q_3(i,ny-3,k); Q_3(i,2,k) = Q_3(i,ny-2,k); Q_3(i,ny-1,k) = Q_3(i,3,k); Q_3(i,ny,k) = Q_3(i,4,k)
        Q_4(i,1,k) = Q_4(i,ny-3,k); Q_4(i,2,k) = Q_4(i,ny-2,k); Q_4(i,ny-1,k) = Q_4(i,3,k); Q_4(i,ny,k) = Q_4(i,4,k)
        Q_5(i,1,k) = Q_5(i,ny-3,k); Q_5(i,2,k) = Q_5(i,ny-2,k); Q_5(i,ny-1,k) = Q_5(i,3,k); Q_5(i,ny,k) = Q_5(i,4,k)
    enddo;enddo
    !$cuf kernel do(1) <<<*,*>>>
    do k = 3, nz-2
      Q_1(1,1,k) = Q_1(nx-3,ny-3,k); Q_1(1,2,k) = Q_1(nx-3,ny-2,k); Q_1(2,1,k) = Q_1(nx-2,ny-3,k); Q_1(2,2,k) = Q_1(nx-2,ny-2,k)
      Q_1(nx-1,1,k) = Q_1(3,ny-3,k); Q_1(nx-1,2,k) = Q_1(3,ny-2,k); Q_1(nx,1,k) = Q_1(4,ny-3,k); Q_1(nx,2,k) = Q_1(4,ny-2,k)
      Q_1(1,ny-1,k) = Q_1(nx-3,3,k); Q_1(1,ny,k) = Q_1(nx-3,4,k); Q_1(2,ny-1,k) = Q_1(nx-2,3,k); Q_1(2,ny,k) = Q_1(nx-2,4,k)
      Q_1(nx-1,ny-1,k) = Q_1(3,3,k); Q_1(nx-1,ny,k) = Q_1(3,4,k); Q_1(nx,ny-1,k) = Q_1(4,3,k); Q_1(nx,ny,k) = Q_1(4,4,k)

      Q_2(1,1,k) = Q_2(nx-3,ny-3,k); Q_2(1,2,k) = Q_2(nx-3,ny-2,k); Q_2(2,1,k) = Q_2(nx-2,ny-3,k); Q_2(2,2,k) = Q_2(nx-2,ny-2,k)
      Q_2(nx-1,1,k) = Q_2(3,ny-3,k); Q_2(nx-1,2,k) = Q_2(3,ny-2,k); Q_2(nx,1,k) = Q_2(4,ny-3,k); Q_2(nx,2,k) = Q_2(4,ny-2,k)
      Q_2(1,ny-1,k) = Q_2(nx-3,3,k); Q_2(1,ny,k) = Q_2(nx-3,4,k); Q_2(2,ny-1,k) = Q_2(nx-2,3,k); Q_2(2,ny,k) = Q_2(nx-2,4,k)
      Q_2(nx-1,ny-1,k) = Q_2(3,3,k); Q_2(nx-1,ny,k) = Q_2(3,4,k); Q_2(nx,ny-1,k) = Q_2(4,3,k); Q_2(nx,ny,k) = Q_2(4,4,k)

      Q_3(1,1,k) = Q_3(nx-3,ny-3,k); Q_3(1,2,k) = Q_3(nx-3,ny-2,k); Q_3(2,1,k) = Q_3(nx-2,ny-3,k); Q_3(2,2,k) = Q_3(nx-2,ny-2,k)
      Q_3(nx-1,1,k) = Q_3(3,ny-3,k); Q_3(nx-1,2,k) = Q_3(3,ny-2,k); Q_3(nx,1,k) = Q_3(4,ny-3,k); Q_3(nx,2,k) = Q_3(4,ny-2,k)
      Q_3(1,ny-1,k) = Q_3(nx-3,3,k); Q_3(1,ny,k) = Q_3(nx-3,4,k); Q_3(2,ny-1,k) = Q_3(nx-2,3,k); Q_3(2,ny,k) = Q_3(nx-2,4,k)
      Q_3(nx-1,ny-1,k) = Q_3(3,3,k); Q_3(nx-1,ny,k) = Q_3(3,4,k); Q_3(nx,ny-1,k) = Q_3(4,3,k); Q_3(nx,ny,k) = Q_3(4,4,k)

      Q_4(1,1,k) = Q_4(nx-3,ny-3,k); Q_4(1,2,k) = Q_4(nx-3,ny-2,k); Q_4(2,1,k) = Q_4(nx-2,ny-3,k); Q_4(2,2,k) = Q_4(nx-2,ny-2,k)
      Q_4(nx-1,1,k) = Q_4(3,ny-3,k); Q_4(nx-1,2,k) = Q_4(3,ny-2,k); Q_4(nx,1,k) = Q_4(4,ny-3,k); Q_4(nx,2,k) = Q_4(4,ny-2,k)
      Q_4(1,ny-1,k) = Q_4(nx-3,3,k); Q_4(1,ny,k) = Q_4(nx-3,4,k); Q_4(2,ny-1,k) = Q_4(nx-2,3,k); Q_4(2,ny,k) = Q_4(nx-2,4,k)
      Q_4(nx-1,ny-1,k) = Q_4(3,3,k); Q_4(nx-1,ny,k) = Q_4(3,4,k); Q_4(nx,ny-1,k) = Q_4(4,3,k); Q_4(nx,ny,k) = Q_4(4,4,k)

      Q_5(1,1,k) = Q_5(nx-3,ny-3,k); Q_5(1,2,k) = Q_5(nx-3,ny-2,k); Q_5(2,1,k) = Q_5(nx-2,ny-3,k); Q_5(2,2,k) = Q_5(nx-2,ny-2,k)
      Q_5(nx-1,1,k) = Q_5(3,ny-3,k); Q_5(nx-1,2,k) = Q_5(3,ny-2,k); Q_5(nx,1,k) = Q_5(4,ny-3,k); Q_5(nx,2,k) = Q_5(4,ny-2,k)
      Q_5(1,ny-1,k) = Q_5(nx-3,3,k); Q_5(1,ny,k) = Q_5(nx-3,4,k); Q_5(2,ny-1,k) = Q_5(nx-2,3,k); Q_5(2,ny,k) = Q_5(nx-2,4,k)
      Q_5(nx-1,ny-1,k) = Q_5(3,3,k); Q_5(nx-1,ny,k) = Q_5(3,4,k); Q_5(nx,ny-1,k) = Q_5(4,3,k); Q_5(nx,ny,k) = Q_5(4,4,k)
    enddo
    !$cuf kernel do(2) <<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        Q_1(i,j,1) = Q_1(i,j,nz-3); Q_1(i,j,2) = Q_1(i,j,nz-2); Q_1(i,j,nz-1) = Q_1(i,j,3); Q_1(i,j,nz) = Q_1(i,j,4)
        Q_2(i,j,1) = Q_2(i,j,nz-3); Q_2(i,j,2) = Q_2(i,j,nz-2); Q_2(i,j,nz-1) = Q_2(i,j,3); Q_2(i,j,nz) = Q_2(i,j,4)
        Q_3(i,j,1) = Q_3(i,j,nz-3); Q_3(i,j,2) = Q_3(i,j,nz-2); Q_3(i,j,nz-1) = Q_3(i,j,3); Q_3(i,j,nz) = Q_3(i,j,4)
        Q_4(i,j,1) = Q_4(i,j,nz-3); Q_4(i,j,2) = Q_4(i,j,nz-2); Q_4(i,j,nz-1) = Q_4(i,j,3); Q_4(i,j,nz) = Q_4(i,j,4)
        Q_5(i,j,1) = Q_5(i,j,nz-3); Q_5(i,j,2) = Q_5(i,j,nz-2); Q_5(i,j,nz-1) = Q_5(i,j,3); Q_5(i,j,nz) = Q_5(i,j,4)
    enddo;enddo
  end subroutine set_bc_cyclic4

  !> Cyclic boundary condition for 6th-order accuracy
  !> Executed on GPU during main time-stepping loop
  subroutine set_bc_cyclic6(id_accuracy, nx, ny, nz, Q_1, Q_2, Q_3, Q_4, Q_5)
    integer(kind=8), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout), device     :: Q_1(nx,ny,nz), Q_2(nx,ny,nz), Q_3(nx,ny,nz), Q_4(nx,ny,nz), Q_5(nx,ny,nz)
    integer i, j, k
    !$cuf kernel do(2)<<<*,*>>>
    do k = 4, nz-3
      do j = 4, ny-3
        Q_1(1,j,k) = Q_1(nx-5,j,k); Q_1(2,j,k) = Q_1(nx-4,j,k); Q_1(3,j,k) = Q_1(nx-3,j,k)
        Q_1(nx-2,j,k) = Q_1(4,j,k); Q_1(nx-1,j,k) = Q_1(5,j,k); Q_1(nx,j,k)   = Q_1(6,j,k)
        Q_2(1,j,k) = Q_2(nx-5,j,k); Q_2(2,j,k) = Q_2(nx-4,j,k); Q_2(3,j,k) = Q_2(nx-3,j,k)
        Q_2(nx-2,j,k) = Q_2(4,j,k); Q_2(nx-1,j,k) = Q_2(5,j,k); Q_2(nx,j,k)   = Q_2(6,j,k)
        Q_3(1,j,k) = Q_3(nx-5,j,k); Q_3(2,j,k) = Q_3(nx-4,j,k); Q_3(3,j,k) = Q_3(nx-3,j,k)
        Q_3(nx-2,j,k) = Q_3(4,j,k); Q_3(nx-1,j,k) = Q_3(5,j,k); Q_3(nx,j,k)   = Q_3(6,j,k)
        Q_4(1,j,k) = Q_4(nx-5,j,k); Q_4(2,j,k) = Q_4(nx-4,j,k); Q_4(3,j,k) = Q_4(nx-3,j,k)
        Q_4(nx-2,j,k) = Q_4(4,j,k); Q_4(nx-1,j,k) = Q_4(5,j,k); Q_4(nx,j,k)   = Q_4(6,j,k)
        Q_5(1,j,k) = Q_5(nx-5,j,k); Q_5(2,j,k) = Q_5(nx-4,j,k); Q_5(3,j,k) = Q_5(nx-3,j,k)
        Q_5(nx-2,j,k) = Q_5(4,j,k); Q_5(nx-1,j,k) = Q_5(5,j,k); Q_5(nx,j,k)   = Q_5(6,j,k)
    enddo;enddo
    !$cuf kernel do(2)<<<*,*>>>
    do k = 4, nz-3
      do i = 4, nx-3
        Q_1(i,1,k) = Q_1(i,ny-5,k); Q_1(i,2,k) = Q_1(i,ny-4,k); Q_1(i,3,k) = Q_1(i,ny-3,k)
        Q_1(i,ny-2,k) = Q_1(i,4,k); Q_1(i,ny-1,k) = Q_1(i,5,k); Q_1(i,ny,k)   = Q_1(i,6,k)
        Q_2(i,1,k) = Q_2(i,ny-5,k); Q_2(i,2,k) = Q_2(i,ny-4,k); Q_2(i,3,k) = Q_2(i,ny-3,k)
        Q_2(i,ny-2,k) = Q_2(i,4,k); Q_2(i,ny-1,k) = Q_2(i,5,k); Q_2(i,ny,k)   = Q_2(i,6,k)
        Q_3(i,1,k) = Q_3(i,ny-5,k); Q_3(i,2,k) = Q_3(i,ny-4,k); Q_3(i,3,k) = Q_3(i,ny-3,k)
        Q_3(i,ny-2,k) = Q_3(i,4,k); Q_3(i,ny-1,k) = Q_3(i,5,k); Q_3(i,ny,k)   = Q_3(i,6,k)
        Q_4(i,1,k) = Q_4(i,ny-5,k); Q_4(i,2,k) = Q_4(i,ny-4,k); Q_4(i,3,k) = Q_4(i,ny-3,k)
        Q_4(i,ny-2,k) = Q_4(i,4,k); Q_4(i,ny-1,k) = Q_4(i,5,k); Q_4(i,ny,k)   = Q_4(i,6,k)
        Q_5(i,1,k) = Q_5(i,ny-5,k); Q_5(i,2,k) = Q_5(i,ny-4,k); Q_5(i,3,k) = Q_5(i,ny-3,k)
        Q_5(i,ny-2,k) = Q_5(i,4,k); Q_5(i,ny-1,k) = Q_5(i,5,k); Q_5(i,ny,k)   = Q_5(i,6,k)
    enddo;enddo
    !$cuf kernel do(1)<<<*,*>>>
    do k = 4, nz-3
      Q_1(1,1,k) = Q_1(nx-5,ny-5,k); Q_1(1,2,k) = Q_1(nx-5,ny-4,k); Q_1(1,3,k) = Q_1(nx-5,ny-3,k)
      Q_1(2,1,k) = Q_1(nx-4,ny-5,k); Q_1(2,2,k) = Q_1(nx-4,ny-4,k); Q_1(2,3,k) = Q_1(nx-4,ny-3,k)
      Q_1(3,1,k) = Q_1(nx-3,ny-5,k); Q_1(3,2,k) = Q_1(nx-3,ny-4,k); Q_1(3,3,k) = Q_1(nx-3,ny-3,k)
      Q_1(nx-2,1,k) = Q_1(4,ny-5,k); Q_1(nx-2,2,k) = Q_1(4,ny-4,k); Q_1(nx-2,3,k) = Q_1(4,ny-3,k)
      Q_1(nx-1,1,k) = Q_1(5,ny-5,k); Q_1(nx-1,2,k) = Q_1(5,ny-4,k); Q_1(nx-1,3,k) = Q_1(5,ny-3,k)
      Q_1(nx,1,k)   = Q_1(6,ny-5,k); Q_1(nx,2,k)   = Q_1(6,ny-4,k); Q_1(nx,3,k)   = Q_1(6,ny-3,k)
      Q_1(1,ny-2,k) = Q_1(nx-5,4,k); Q_1(1,ny-1,k) = Q_1(nx-5,5,k); Q_1(1,ny,k)   = Q_1(nx-5,6,k)
      Q_1(2,ny-2,k) = Q_1(nx-4,4,k); Q_1(2,ny-1,k) = Q_1(nx-4,5,k); Q_1(2,ny,k)   = Q_1(nx-4,6,k)
      Q_1(3,ny-2,k) = Q_1(nx-3,4,k); Q_1(3,ny-1,k) = Q_1(nx-3,5,k); Q_1(3,ny,k)   = Q_1(nx-3,6,k)
      Q_1(nx-2,ny-2,k) = Q_1(4,4,k); Q_1(nx-2,ny-1,k) = Q_1(4,5,k); Q_1(nx-2,ny,k)   = Q_1(4,6,k)
      Q_1(nx-1,ny-2,k) = Q_1(5,4,k); Q_1(nx-1,ny-1,k) = Q_1(5,5,k); Q_1(nx-1,ny,k)   = Q_1(5,6,k)
      Q_1(nx,ny-2,k)   = Q_1(6,4,k); Q_1(nx,ny-1,k)   = Q_1(6,5,k); Q_1(nx,ny,k)     = Q_1(6,6,k)

      Q_2(1,1,k) = Q_2(nx-5,ny-5,k); Q_2(1,2,k) = Q_2(nx-5,ny-4,k); Q_2(1,3,k) = Q_2(nx-5,ny-3,k)
      Q_2(2,1,k) = Q_2(nx-4,ny-5,k); Q_2(2,2,k) = Q_2(nx-4,ny-4,k); Q_2(2,3,k) = Q_2(nx-4,ny-3,k)
      Q_2(3,1,k) = Q_2(nx-3,ny-5,k); Q_2(3,2,k) = Q_2(nx-3,ny-4,k); Q_2(3,3,k) = Q_2(nx-3,ny-3,k)
      Q_2(nx-2,1,k) = Q_2(4,ny-5,k); Q_2(nx-2,2,k) = Q_2(4,ny-4,k); Q_2(nx-2,3,k) = Q_2(4,ny-3,k)
      Q_2(nx-1,1,k) = Q_2(5,ny-5,k); Q_2(nx-1,2,k) = Q_2(5,ny-4,k); Q_2(nx-1,3,k) = Q_2(5,ny-3,k)
      Q_2(nx,1,k)   = Q_2(6,ny-5,k); Q_2(nx,2,k)   = Q_2(6,ny-4,k); Q_2(nx,3,k)   = Q_2(6,ny-3,k)
      Q_2(1,ny-2,k) = Q_2(nx-5,4,k); Q_2(1,ny-1,k) = Q_2(nx-5,5,k); Q_2(1,ny,k)   = Q_2(nx-5,6,k)
      Q_2(2,ny-2,k) = Q_2(nx-4,4,k); Q_2(2,ny-1,k) = Q_2(nx-4,5,k); Q_2(2,ny,k)   = Q_2(nx-4,6,k)
      Q_2(3,ny-2,k) = Q_2(nx-3,4,k); Q_2(3,ny-1,k) = Q_2(nx-3,5,k); Q_2(3,ny,k)   = Q_2(nx-3,6,k)
      Q_2(nx-2,ny-2,k) = Q_2(4,4,k); Q_2(nx-2,ny-1,k) = Q_2(4,5,k); Q_2(nx-2,ny,k)   = Q_2(4,6,k)
      Q_2(nx-1,ny-2,k) = Q_2(5,4,k); Q_2(nx-1,ny-1,k) = Q_2(5,5,k); Q_2(nx-1,ny,k)   = Q_2(5,6,k)
      Q_2(nx,ny-2,k)   = Q_2(6,4,k); Q_2(nx,ny-1,k)   = Q_2(6,5,k); Q_2(nx,ny,k)     = Q_2(6,6,k)

      Q_3(1,1,k) = Q_3(nx-5,ny-5,k); Q_3(1,2,k) = Q_3(nx-5,ny-4,k); Q_3(1,3,k) = Q_3(nx-5,ny-3,k)
      Q_3(2,1,k) = Q_3(nx-4,ny-5,k); Q_3(2,2,k) = Q_3(nx-4,ny-4,k); Q_3(2,3,k) = Q_3(nx-4,ny-3,k)
      Q_3(3,1,k) = Q_3(nx-3,ny-5,k); Q_3(3,2,k) = Q_3(nx-3,ny-4,k); Q_3(3,3,k) = Q_3(nx-3,ny-3,k)
      Q_3(nx-2,1,k) = Q_3(4,ny-5,k); Q_3(nx-2,2,k) = Q_3(4,ny-4,k); Q_3(nx-2,3,k) = Q_3(4,ny-3,k)
      Q_3(nx-1,1,k) = Q_3(5,ny-5,k); Q_3(nx-1,2,k) = Q_3(5,ny-4,k); Q_3(nx-1,3,k) = Q_3(5,ny-3,k)
      Q_3(nx,1,k)   = Q_3(6,ny-5,k); Q_3(nx,2,k)   = Q_3(6,ny-4,k); Q_3(nx,3,k)   = Q_3(6,ny-3,k)
      Q_3(1,ny-2,k) = Q_3(nx-5,4,k); Q_3(1,ny-1,k) = Q_3(nx-5,5,k); Q_3(1,ny,k)   = Q_3(nx-5,6,k)
      Q_3(2,ny-2,k) = Q_3(nx-4,4,k); Q_3(2,ny-1,k) = Q_3(nx-4,5,k); Q_3(2,ny,k)   = Q_3(nx-4,6,k)
      Q_3(3,ny-2,k) = Q_3(nx-3,4,k); Q_3(3,ny-1,k) = Q_3(nx-3,5,k); Q_3(3,ny,k)   = Q_3(nx-3,6,k)
      Q_3(nx-2,ny-2,k) = Q_3(4,4,k); Q_3(nx-2,ny-1,k) = Q_3(4,5,k); Q_3(nx-2,ny,k)   = Q_3(4,6,k)
      Q_3(nx-1,ny-2,k) = Q_3(5,4,k); Q_3(nx-1,ny-1,k) = Q_3(5,5,k); Q_3(nx-1,ny,k)   = Q_3(5,6,k)
      Q_3(nx,ny-2,k)   = Q_3(6,4,k); Q_3(nx,ny-1,k)   = Q_3(6,5,k); Q_3(nx,ny,k)     = Q_3(6,6,k)

      Q_4(1,1,k) = Q_4(nx-5,ny-5,k); Q_4(1,2,k) = Q_4(nx-5,ny-4,k); Q_4(1,3,k) = Q_4(nx-5,ny-3,k)
      Q_4(2,1,k) = Q_4(nx-4,ny-5,k); Q_4(2,2,k) = Q_4(nx-4,ny-4,k); Q_4(2,3,k) = Q_4(nx-4,ny-3,k)
      Q_4(3,1,k) = Q_4(nx-3,ny-5,k); Q_4(3,2,k) = Q_4(nx-3,ny-4,k); Q_4(3,3,k) = Q_4(nx-3,ny-3,k)
      Q_4(nx-2,1,k) = Q_4(4,ny-5,k); Q_4(nx-2,2,k) = Q_4(4,ny-4,k); Q_4(nx-2,3,k) = Q_4(4,ny-3,k)
      Q_4(nx-1,1,k) = Q_4(5,ny-5,k); Q_4(nx-1,2,k) = Q_4(5,ny-4,k); Q_4(nx-1,3,k) = Q_4(5,ny-3,k)
      Q_4(nx,1,k)   = Q_4(6,ny-5,k); Q_4(nx,2,k)   = Q_4(6,ny-4,k); Q_4(nx,3,k)   = Q_4(6,ny-3,k)
      Q_4(1,ny-2,k) = Q_4(nx-5,4,k); Q_4(1,ny-1,k) = Q_4(nx-5,5,k); Q_4(1,ny,k)   = Q_4(nx-5,6,k)
      Q_4(2,ny-2,k) = Q_4(nx-4,4,k); Q_4(2,ny-1,k) = Q_4(nx-4,5,k); Q_4(2,ny,k)   = Q_4(nx-4,6,k)
      Q_4(3,ny-2,k) = Q_4(nx-3,4,k); Q_4(3,ny-1,k) = Q_4(nx-3,5,k); Q_4(3,ny,k)   = Q_4(nx-3,6,k)
      Q_4(nx-2,ny-2,k) = Q_4(4,4,k); Q_4(nx-2,ny-1,k) = Q_4(4,5,k); Q_4(nx-2,ny,k)   = Q_4(4,6,k)
      Q_4(nx-1,ny-2,k) = Q_4(5,4,k); Q_4(nx-1,ny-1,k) = Q_4(5,5,k); Q_4(nx-1,ny,k)   = Q_4(5,6,k)
      Q_4(nx,ny-2,k)   = Q_4(6,4,k); Q_4(nx,ny-1,k)   = Q_4(6,5,k); Q_4(nx,ny,k)     = Q_4(6,6,k)

      Q_5(1,1,k) = Q_5(nx-5,ny-5,k); Q_5(1,2,k) = Q_5(nx-5,ny-4,k); Q_5(1,3,k) = Q_5(nx-5,ny-3,k)
      Q_5(2,1,k) = Q_5(nx-4,ny-5,k); Q_5(2,2,k) = Q_5(nx-4,ny-4,k); Q_5(2,3,k) = Q_5(nx-4,ny-3,k)
      Q_5(3,1,k) = Q_5(nx-3,ny-5,k); Q_5(3,2,k) = Q_5(nx-3,ny-4,k); Q_5(3,3,k) = Q_5(nx-3,ny-3,k)
      Q_5(nx-2,1,k) = Q_5(4,ny-5,k); Q_5(nx-2,2,k) = Q_5(4,ny-4,k); Q_5(nx-2,3,k) = Q_5(4,ny-3,k)
      Q_5(nx-1,1,k) = Q_5(5,ny-5,k); Q_5(nx-1,2,k) = Q_5(5,ny-4,k); Q_5(nx-1,3,k) = Q_5(5,ny-3,k)
      Q_5(nx,1,k)   = Q_5(6,ny-5,k); Q_5(nx,2,k)   = Q_5(6,ny-4,k); Q_5(nx,3,k)   = Q_5(6,ny-3,k)
      Q_5(1,ny-2,k) = Q_5(nx-5,4,k); Q_5(1,ny-1,k) = Q_5(nx-5,5,k); Q_5(1,ny,k)   = Q_5(nx-5,6,k)
      Q_5(2,ny-2,k) = Q_5(nx-4,4,k); Q_5(2,ny-1,k) = Q_5(nx-4,5,k); Q_5(2,ny,k)   = Q_5(nx-4,6,k)
      Q_5(3,ny-2,k) = Q_5(nx-3,4,k); Q_5(3,ny-1,k) = Q_5(nx-3,5,k); Q_5(3,ny,k)   = Q_5(nx-3,6,k)
      Q_5(nx-2,ny-2,k) = Q_5(4,4,k); Q_5(nx-2,ny-1,k) = Q_5(4,5,k); Q_5(nx-2,ny,k)   = Q_5(4,6,k)
      Q_5(nx-1,ny-2,k) = Q_5(5,4,k); Q_5(nx-1,ny-1,k) = Q_5(5,5,k); Q_5(nx-1,ny,k)   = Q_5(5,6,k)
      Q_5(nx,ny-2,k)   = Q_5(6,4,k); Q_5(nx,ny-1,k)   = Q_5(6,5,k); Q_5(nx,ny,k)     = Q_5(6,6,k)
    enddo
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        Q_1(i,j,1) = Q_1(i,j,nz-5); Q_1(i,j,2) = Q_1(i,j,nz-4); Q_1(i,j,3) = Q_1(i,j,nz-3)
        Q_1(i,j,nz-2) = Q_1(i,j,4); Q_1(i,j,nz-1) = Q_1(i,j,5); Q_1(i,j,nz)   = Q_1(i,j,6)
        Q_2(i,j,1) = Q_2(i,j,nz-5); Q_2(i,j,2) = Q_2(i,j,nz-4); Q_2(i,j,3) = Q_2(i,j,nz-3)
        Q_2(i,j,nz-2) = Q_2(i,j,4); Q_2(i,j,nz-1) = Q_2(i,j,5); Q_2(i,j,nz)   = Q_2(i,j,6)
        Q_3(i,j,1) = Q_3(i,j,nz-5); Q_3(i,j,2) = Q_3(i,j,nz-4); Q_3(i,j,3) = Q_3(i,j,nz-3)
        Q_3(i,j,nz-2) = Q_3(i,j,4); Q_3(i,j,nz-1) = Q_3(i,j,5); Q_3(i,j,nz)   = Q_3(i,j,6)
        Q_4(i,j,1) = Q_4(i,j,nz-5); Q_4(i,j,2) = Q_4(i,j,nz-4); Q_4(i,j,3) = Q_4(i,j,nz-3)
        Q_4(i,j,nz-2) = Q_4(i,j,4); Q_4(i,j,nz-1) = Q_4(i,j,5); Q_4(i,j,nz)   = Q_4(i,j,6)
        Q_5(i,j,1) = Q_5(i,j,nz-5); Q_5(i,j,2) = Q_5(i,j,nz-4); Q_5(i,j,3) = Q_5(i,j,nz-3)
        Q_5(i,j,nz-2) = Q_5(i,j,4); Q_5(i,j,nz-1) = Q_5(i,j,5); Q_5(i,j,nz)   = Q_5(i,j,6)
    enddo;enddo
  end subroutine set_bc_cyclic6

  !> Cyclic boundary condition in z-direction for 6th-order accuracy
  !> Executed on GPU during main time-stepping loop
  subroutine set_bc_cyclic_z(nx, ny, nz, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(inout), device :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    integer i, j
    !$cuf kernel do(2)<<<*,*>>>
    do j = 1, ny
      do i = 1, nx
        QJ_1(i,j,1) = QJ_1(i,j,nz-5); QJ_1(i,j,2) = QJ_1(i,j,nz-4); QJ_1(i,j,3) = QJ_1(i,j,nz-3)
        QJ_1(i,j,nz-2) = QJ_1(i,j,4); QJ_1(i,j,nz-1) = QJ_1(i,j,5); QJ_1(i,j,nz)   = QJ_1(i,j,6)
        QJ_2(i,j,1) = QJ_2(i,j,nz-5); QJ_2(i,j,2) = QJ_2(i,j,nz-4); QJ_2(i,j,3) = QJ_2(i,j,nz-3)
        QJ_2(i,j,nz-2) = QJ_2(i,j,4); QJ_2(i,j,nz-1) = QJ_2(i,j,5); QJ_2(i,j,nz)   = QJ_2(i,j,6)
        QJ_3(i,j,1) = QJ_3(i,j,nz-5); QJ_3(i,j,2) = QJ_3(i,j,nz-4); QJ_3(i,j,3) = QJ_3(i,j,nz-3)
        QJ_3(i,j,nz-2) = QJ_3(i,j,4); QJ_3(i,j,nz-1) = QJ_3(i,j,5); QJ_3(i,j,nz)   = QJ_3(i,j,6)
        QJ_4(i,j,1) = QJ_4(i,j,nz-5); QJ_4(i,j,2) = QJ_4(i,j,nz-4); QJ_4(i,j,3) = QJ_4(i,j,nz-3)
        QJ_4(i,j,nz-2) = QJ_4(i,j,4); QJ_4(i,j,nz-1) = QJ_4(i,j,5); QJ_4(i,j,nz)   = QJ_4(i,j,6)
        QJ_5(i,j,1) = QJ_5(i,j,nz-5); QJ_5(i,j,2) = QJ_5(i,j,nz-4); QJ_5(i,j,3) = QJ_5(i,j,nz-3)
        QJ_5(i,j,nz-2) = QJ_5(i,j,4); QJ_5(i,j,nz-1) = QJ_5(i,j,5); QJ_5(i,j,nz)   = QJ_5(i,j,6)
    enddo;enddo
  end subroutine set_bc_cyclic_z

  !> Boundary condition for SGS viscosity and turbulent kinetic energy
  !> Executed on GPU during main time-stepping loop
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

