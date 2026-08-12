module set_bc_common
  use mod_globals, only : nx, ny
  implicit none
  interface set_bc_cyclic
    module procedure set_bc_cyclic2_init, set_bc_cyclic2, set_bc_cyclic4, &
                     set_bc_cyclic4_init, set_bc_cyclic6, set_bc_cyclic6_init
  end interface set_bc_cyclic
contains
  !> Cyclic boundary condition initialization for second-order accuracy
  !> Executed on CPU before main time-stepping loop
  subroutine set_bc_cyclic2_init(id_accuracy, nx, ny, Q)
    integer(kind=2), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny
    real(8), intent(inout)             :: Q(nx,ny,4)
    integer i, j
    do j = 2, ny-1
      Q(1,j,:)  = Q(nx-1,j,:)
      Q(nx,j,:) = Q(2,j,:)
    enddo
    do i = 2, nx-1
      Q(i,1,:)  = Q(i,ny-1,:)
      Q(i,ny,:) = Q(i,2,:)
    enddo
    Q(1,1,:)   = Q(nx-1,ny-1,:)
    Q(nx,1,:)  = Q(2,ny-1,:)
    Q(1,ny,:)  = Q(nx-1,2,:)
    Q(nx,ny,:) = Q(2,2,:)
  end subroutine set_bc_cyclic2_init

  !> Cyclic boundary condition initialization for 4th-order accuracy
  !> Executed on CPU before main time-stepping loop
  subroutine set_bc_cyclic4_init(id_accuracy, nx, ny, Q)
    integer(kind=4), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny
    real(8), intent(inout)             :: Q(nx,ny,4)
    integer i, j
    do j = 3, ny-2
      Q(1:2,j,:) = Q(nx-3:nx-2,j,:)
      Q(nx-1:nx,j,:) = Q(3:4,j,:)
    enddo
    do i = 3, nx-2
      Q(i,1:2,:) = Q(i,ny-3:ny-2,:)
      Q(i,ny-1:ny,:) = Q(i,3:4,:)
    enddo
    Q(1:2,1:2,:) = Q(nx-3:nx-2,ny-3:ny-2,:)
    Q(nx-1:nx,1:2,:) = Q(3:4,ny-3:ny-2,:)
    Q(1:2,ny-1:ny,:) = Q(nx-3:nx-2,3:4,:)
    Q(nx-1:nx,ny-1:ny,:) = Q(3:4,3:4,:)
  end subroutine set_bc_cyclic4_init

  !> Cyclic boundary condition initialization for 6th-order accuracy
  !> Executed on CPU before main time-stepping loop
  subroutine set_bc_cyclic6_init(id_accuracy, nx, ny, Q)
    integer(kind=8), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny
    real(8), intent(inout)             :: Q(nx,ny,4)
    integer i, j
    do j = 4, ny-3
      Q(1:3,j,:) = Q(nx-5:nx-3,j,:)
      Q(nx-2:nx,j,:) = Q(4:6,j,:)
    enddo
    do i = 4, nx-3
      Q(i,1:3,:) = Q(i,ny-5:ny-3,:)
      Q(i,ny-2:ny,:) = Q(i,4:6,:)
    enddo
    Q(1:3,1:3,:) = Q(nx-5:nx-3,ny-5:ny-3,:)
    Q(nx-2:nx,1:3,:) = Q(4:6,ny-5:ny-3,:)
    Q(1:3,ny-2:ny,:) = Q(nx-5:nx-3,4:6,:)
    Q(nx-2:nx,ny-2:ny,:) = Q(4:6,4:6,:)
  end subroutine set_bc_cyclic6_init

  !> Cyclic boundary condition for second-order accuracy (SoA, one component per call)
  !> Executed on GPU during main time-stepping loop
  subroutine set_bc_cyclic2(id_accuracy, nx, ny, Q_1, Q_2, Q_3, Q_4)
    integer(kind=2), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny
    real(8), intent(inout), device     :: Q_1(nx,ny), Q_2(nx,ny), Q_3(nx,ny), Q_4(nx,ny)
    integer i, j
    !$cuf kernel do(1)<<<*,*>>>
    do j = 2, ny-1
      Q_1(1,j)  = Q_1(nx-1,j); Q_2(1,j)  = Q_2(nx-1,j); Q_3(1,j)  = Q_3(nx-1,j); Q_4(1,j)  = Q_4(nx-1,j)
      Q_1(nx,j) = Q_1(2,j);    Q_2(nx,j) = Q_2(2,j);    Q_3(nx,j) = Q_3(2,j);    Q_4(nx,j) = Q_4(2,j)
    enddo
    !$cuf kernel do(1)<<<*,*>>>
    do i = 2, nx-1
      Q_1(i,1)  = Q_1(i,ny-1); Q_2(i,1)  = Q_2(i,ny-1); Q_3(i,1)  = Q_3(i,ny-1); Q_4(i,1)  = Q_4(i,ny-1)
      Q_1(i,ny) = Q_1(i,2);    Q_2(i,ny) = Q_2(i,2);    Q_3(i,ny) = Q_3(i,2);    Q_4(i,ny) = Q_4(i,2)
    enddo
    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, 1
      Q_1(1,1)   = Q_1(nx-1,ny-1); Q_2(1,1)   = Q_2(nx-1,ny-1); Q_3(1,1)   = Q_3(nx-1,ny-1); Q_4(1,1)   = Q_4(nx-1,ny-1)
      Q_1(nx,1)  = Q_1(2,ny-1);    Q_2(nx,1)  = Q_2(2,ny-1);    Q_3(nx,1)  = Q_3(2,ny-1);    Q_4(nx,1)  = Q_4(2,ny-1)
      Q_1(1,ny)  = Q_1(nx-1,2);    Q_2(1,ny)  = Q_2(nx-1,2);    Q_3(1,ny)  = Q_3(nx-1,2);    Q_4(1,ny)  = Q_4(nx-1,2)
      Q_1(nx,ny) = Q_1(2,2);       Q_2(nx,ny) = Q_2(2,2);       Q_3(nx,ny) = Q_3(2,2);       Q_4(nx,ny) = Q_4(2,2)
    enddo
  end subroutine set_bc_cyclic2

  !> Cyclic boundary condition for 4th-order accuracy (SoA, one component per call)
  !> Executed on GPU during main time-stepping loop
  subroutine set_bc_cyclic4(id_accuracy, nx, ny, Q_1, Q_2, Q_3, Q_4)
    integer(kind=4), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny
    real(8), intent(inout), device     :: Q_1(nx,ny), Q_2(nx,ny), Q_3(nx,ny), Q_4(nx,ny)
    integer i, j
    !$cuf kernel do(1) <<<*,*>>>
    do j = 3, ny-2
      Q_1(1,j) = Q_1(nx-3,j); Q_2(1,j) = Q_2(nx-3,j); Q_3(1,j) = Q_3(nx-3,j); Q_4(1,j) = Q_4(nx-3,j)
      Q_1(2,j) = Q_1(nx-2,j); Q_2(2,j) = Q_2(nx-2,j); Q_3(2,j) = Q_3(nx-2,j); Q_4(2,j) = Q_4(nx-2,j)
      Q_1(nx-1,j) = Q_1(3,j); Q_2(nx-1,j) = Q_2(3,j); Q_3(nx-1,j) = Q_3(3,j); Q_4(nx-1,j) = Q_4(3,j)
      Q_1(nx,j)   = Q_1(4,j); Q_2(nx,j)   = Q_2(4,j); Q_3(nx,j)   = Q_3(4,j); Q_4(nx,j)   = Q_4(4,j)
    enddo
    !$cuf kernel do(1) <<<*,*>>>
    do i = 3, nx-2
      Q_1(i,1) = Q_1(i,ny-3); Q_2(i,1) = Q_2(i,ny-3); Q_3(i,1) = Q_3(i,ny-3); Q_4(i,1) = Q_4(i,ny-3)
      Q_1(i,2) = Q_1(i,ny-2); Q_2(i,2) = Q_2(i,ny-2); Q_3(i,2) = Q_3(i,ny-2); Q_4(i,2) = Q_4(i,ny-2)
      Q_1(i,ny-1) = Q_1(i,3); Q_2(i,ny-1) = Q_2(i,3); Q_3(i,ny-1) = Q_3(i,3); Q_4(i,ny-1) = Q_4(i,3)
      Q_1(i,ny)   = Q_1(i,4); Q_2(i,ny)   = Q_2(i,4); Q_3(i,ny)   = Q_3(i,4); Q_4(i,ny)   = Q_4(i,4)
    enddo
    !$cuf kernel do(1) <<<*,*>>>
    do i = 1, 1
      Q_1(1,1) = Q_1(nx-3,ny-3); Q_2(1,1) = Q_2(nx-3,ny-3); Q_3(1,1) = Q_3(nx-3,ny-3); Q_4(1,1) = Q_4(nx-3,ny-3)
      Q_1(1,2) = Q_1(nx-3,ny-2); Q_2(1,2) = Q_2(nx-3,ny-2); Q_3(1,2) = Q_3(nx-3,ny-2); Q_4(1,2) = Q_4(nx-3,ny-2)
      Q_1(2,1) = Q_1(nx-2,ny-3); Q_2(2,1) = Q_2(nx-2,ny-3); Q_3(2,1) = Q_3(nx-2,ny-3); Q_4(2,1) = Q_4(nx-2,ny-3)
      Q_1(2,2) = Q_1(nx-2,ny-2); Q_2(2,2) = Q_2(nx-2,ny-2); Q_3(2,2) = Q_3(nx-2,ny-2); Q_4(2,2) = Q_4(nx-2,ny-2)
      Q_1(nx-1,1) = Q_1(3,ny-3); Q_2(nx-1,1) = Q_2(3,ny-3); Q_3(nx-1,1) = Q_3(3,ny-3); Q_4(nx-1,1) = Q_4(3,ny-3)
      Q_1(nx-1,2) = Q_1(3,ny-2); Q_2(nx-1,2) = Q_2(3,ny-2); Q_3(nx-1,2) = Q_3(3,ny-2); Q_4(nx-1,2) = Q_4(3,ny-2)
      Q_1(nx,1)   = Q_1(4,ny-3); Q_2(nx,1)   = Q_2(4,ny-3); Q_3(nx,1)   = Q_3(4,ny-3); Q_4(nx,1)   = Q_4(4,ny-3)
      Q_1(nx,2)   = Q_1(4,ny-2); Q_2(nx,2)   = Q_2(4,ny-2); Q_3(nx,2)   = Q_3(4,ny-2); Q_4(nx,2)   = Q_4(4,ny-2)
      Q_1(1,ny-1) = Q_1(nx-3,3); Q_2(1,ny-1) = Q_2(nx-3,3); Q_3(1,ny-1) = Q_3(nx-3,3); Q_4(1,ny-1) = Q_4(nx-3,3)
      Q_1(1,ny)   = Q_1(nx-3,4); Q_2(1,ny)   = Q_2(nx-3,4); Q_3(1,ny)   = Q_3(nx-3,4); Q_4(1,ny)   = Q_4(nx-3,4)
      Q_1(2,ny-1) = Q_1(nx-2,3); Q_2(2,ny-1) = Q_2(nx-2,3); Q_3(2,ny-1) = Q_3(nx-2,3); Q_4(2,ny-1) = Q_4(nx-2,3)
      Q_1(2,ny)   = Q_1(nx-2,4); Q_2(2,ny)   = Q_2(nx-2,4); Q_3(2,ny)   = Q_3(nx-2,4); Q_4(2,ny)   = Q_4(nx-2,4)
      Q_1(nx-1,ny-1) = Q_1(3,3); Q_2(nx-1,ny-1) = Q_2(3,3); Q_3(nx-1,ny-1) = Q_3(3,3); Q_4(nx-1,ny-1) = Q_4(3,3)
      Q_1(nx-1,ny)   = Q_1(3,4); Q_2(nx-1,ny)   = Q_2(3,4); Q_3(nx-1,ny)   = Q_3(3,4); Q_4(nx-1,ny)   = Q_4(3,4)
      Q_1(nx,ny-1)   = Q_1(4,3); Q_2(nx,ny-1)   = Q_2(4,3); Q_3(nx,ny-1)   = Q_3(4,3); Q_4(nx,ny-1)   = Q_4(4,3)
      Q_1(nx,ny)     = Q_1(4,4); Q_2(nx,ny)     = Q_2(4,4); Q_3(nx,ny)     = Q_3(4,4); Q_4(nx,ny)     = Q_4(4,4)
    enddo
  end subroutine set_bc_cyclic4

  !> Cyclic boundary condition for 6th-order accuracy (SoA, one component per call)
  !> Executed on GPU during main time-stepping loop
  subroutine set_bc_cyclic6(id_accuracy, nx, ny, Q_1, Q_2, Q_3, Q_4)
    integer(kind=8), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny
    real(8), intent(inout), device     :: Q_1(nx,ny), Q_2(nx,ny), Q_3(nx,ny), Q_4(nx,ny)
    integer i, j
    !$cuf kernel do(1)<<<*,*>>>
    do j = 4, ny-3
      Q_1(1,j) = Q_1(nx-5,j); Q_2(1,j) = Q_2(nx-5,j); Q_3(1,j) = Q_3(nx-5,j); Q_4(1,j) = Q_4(nx-5,j)
      Q_1(2,j) = Q_1(nx-4,j); Q_2(2,j) = Q_2(nx-4,j); Q_3(2,j) = Q_3(nx-4,j); Q_4(2,j) = Q_4(nx-4,j)
      Q_1(3,j) = Q_1(nx-3,j); Q_2(3,j) = Q_2(nx-3,j); Q_3(3,j) = Q_3(nx-3,j); Q_4(3,j) = Q_4(nx-3,j)
      Q_1(nx-2,j) = Q_1(4,j); Q_2(nx-2,j) = Q_2(4,j); Q_3(nx-2,j) = Q_3(4,j); Q_4(nx-2,j) = Q_4(4,j)
      Q_1(nx-1,j) = Q_1(5,j); Q_2(nx-1,j) = Q_2(5,j); Q_3(nx-1,j) = Q_3(5,j); Q_4(nx-1,j) = Q_4(5,j)
      Q_1(nx,j)   = Q_1(6,j); Q_2(nx,j)   = Q_2(6,j); Q_3(nx,j)   = Q_3(6,j); Q_4(nx,j)   = Q_4(6,j)
    enddo
    !$cuf kernel do(1)<<<*,*>>>
    do i = 4, nx-3
      Q_1(i,1) = Q_1(i,ny-5); Q_2(i,1) = Q_2(i,ny-5); Q_3(i,1) = Q_3(i,ny-5); Q_4(i,1) = Q_4(i,ny-5)
      Q_1(i,2) = Q_1(i,ny-4); Q_2(i,2) = Q_2(i,ny-4); Q_3(i,2) = Q_3(i,ny-4); Q_4(i,2) = Q_4(i,ny-4)
      Q_1(i,3) = Q_1(i,ny-3); Q_2(i,3) = Q_2(i,ny-3); Q_3(i,3) = Q_3(i,ny-3); Q_4(i,3) = Q_4(i,ny-3)
      Q_1(i,ny-2) = Q_1(i,4); Q_2(i,ny-2) = Q_2(i,4); Q_3(i,ny-2) = Q_3(i,4); Q_4(i,ny-2) = Q_4(i,4)
      Q_1(i,ny-1) = Q_1(i,5); Q_2(i,ny-1) = Q_2(i,5); Q_3(i,ny-1) = Q_3(i,5); Q_4(i,ny-1) = Q_4(i,5)
      Q_1(i,ny)   = Q_1(i,6); Q_2(i,ny)   = Q_2(i,6); Q_3(i,ny)   = Q_3(i,6); Q_4(i,ny)   = Q_4(i,6)
    enddo
    !$cuf kernel do(1)<<<*,*>>>
    do i = 1, 1
      Q_1(1,1) = Q_1(nx-5,ny-5); Q_2(1,1) = Q_2(nx-5,ny-5); Q_3(1,1) = Q_3(nx-5,ny-5); Q_4(1,1) = Q_4(nx-5,ny-5)
      Q_1(1,2) = Q_1(nx-5,ny-4); Q_2(1,2) = Q_2(nx-5,ny-4); Q_3(1,2) = Q_3(nx-5,ny-4); Q_4(1,2) = Q_4(nx-5,ny-4)
      Q_1(1,3) = Q_1(nx-5,ny-3); Q_2(1,3) = Q_2(nx-5,ny-3); Q_3(1,3) = Q_3(nx-5,ny-3); Q_4(1,3) = Q_4(nx-5,ny-3)
      Q_1(2,1) = Q_1(nx-4,ny-5); Q_2(2,1) = Q_2(nx-4,ny-5); Q_3(2,1) = Q_3(nx-4,ny-5); Q_4(2,1) = Q_4(nx-4,ny-5)
      Q_1(2,2) = Q_1(nx-4,ny-4); Q_2(2,2) = Q_2(nx-4,ny-4); Q_3(2,2) = Q_3(nx-4,ny-4); Q_4(2,2) = Q_4(nx-4,ny-4)
      Q_1(2,3) = Q_1(nx-4,ny-3); Q_2(2,3) = Q_2(nx-4,ny-3); Q_3(2,3) = Q_3(nx-4,ny-3); Q_4(2,3) = Q_4(nx-4,ny-3)
      Q_1(3,1) = Q_1(nx-3,ny-5); Q_2(3,1) = Q_2(nx-3,ny-5); Q_3(3,1) = Q_3(nx-3,ny-5); Q_4(3,1) = Q_4(nx-3,ny-5)
      Q_1(3,2) = Q_1(nx-3,ny-4); Q_2(3,2) = Q_2(nx-3,ny-4); Q_3(3,2) = Q_3(nx-3,ny-4); Q_4(3,2) = Q_4(nx-3,ny-4)
      Q_1(3,3) = Q_1(nx-3,ny-3); Q_2(3,3) = Q_2(nx-3,ny-3); Q_3(3,3) = Q_3(nx-3,ny-3); Q_4(3,3) = Q_4(nx-3,ny-3)
      Q_1(nx-2,1) = Q_1(4,ny-5); Q_2(nx-2,1) = Q_2(4,ny-5); Q_3(nx-2,1) = Q_3(4,ny-5); Q_4(nx-2,1) = Q_4(4,ny-5)
      Q_1(nx-2,2) = Q_1(4,ny-4); Q_2(nx-2,2) = Q_2(4,ny-4); Q_3(nx-2,2) = Q_3(4,ny-4); Q_4(nx-2,2) = Q_4(4,ny-4)
      Q_1(nx-2,3) = Q_1(4,ny-3); Q_2(nx-2,3) = Q_2(4,ny-3); Q_3(nx-2,3) = Q_3(4,ny-3); Q_4(nx-2,3) = Q_4(4,ny-3)
      Q_1(nx-1,1) = Q_1(5,ny-5); Q_2(nx-1,1) = Q_2(5,ny-5); Q_3(nx-1,1) = Q_3(5,ny-5); Q_4(nx-1,1) = Q_4(5,ny-5)
      Q_1(nx-1,2) = Q_1(5,ny-4); Q_2(nx-1,2) = Q_2(5,ny-4); Q_3(nx-1,2) = Q_3(5,ny-4); Q_4(nx-1,2) = Q_4(5,ny-4)
      Q_1(nx-1,3) = Q_1(5,ny-3); Q_2(nx-1,3) = Q_2(5,ny-3); Q_3(nx-1,3) = Q_3(5,ny-3); Q_4(nx-1,3) = Q_4(5,ny-3)
      Q_1(nx,1)   = Q_1(6,ny-5); Q_2(nx,1)   = Q_2(6,ny-5); Q_3(nx,1)   = Q_3(6,ny-5); Q_4(nx,1)   = Q_4(6,ny-5)
      Q_1(nx,2)   = Q_1(6,ny-4); Q_2(nx,2)   = Q_2(6,ny-4); Q_3(nx,2)   = Q_3(6,ny-4); Q_4(nx,2)   = Q_4(6,ny-4)
      Q_1(nx,3)   = Q_1(6,ny-3); Q_2(nx,3)   = Q_2(6,ny-3); Q_3(nx,3)   = Q_3(6,ny-3); Q_4(nx,3)   = Q_4(6,ny-3)
      Q_1(1,ny-2) = Q_1(nx-5,4); Q_2(1,ny-2) = Q_2(nx-5,4); Q_3(1,ny-2) = Q_3(nx-5,4); Q_4(1,ny-2) = Q_4(nx-5,4)
      Q_1(1,ny-1) = Q_1(nx-5,5); Q_2(1,ny-1) = Q_2(nx-5,5); Q_3(1,ny-1) = Q_3(nx-5,5); Q_4(1,ny-1) = Q_4(nx-5,5)
      Q_1(1,ny)   = Q_1(nx-5,6); Q_2(1,ny)   = Q_2(nx-5,6); Q_3(1,ny)   = Q_3(nx-5,6); Q_4(1,ny)   = Q_4(nx-5,6)
      Q_1(2,ny-2) = Q_1(nx-4,4); Q_2(2,ny-2) = Q_2(nx-4,4); Q_3(2,ny-2) = Q_3(nx-4,4); Q_4(2,ny-2) = Q_4(nx-4,4)
      Q_1(2,ny-1) = Q_1(nx-4,5); Q_2(2,ny-1) = Q_2(nx-4,5); Q_3(2,ny-1) = Q_3(nx-4,5); Q_4(2,ny-1) = Q_4(nx-4,5)
      Q_1(2,ny)   = Q_1(nx-4,6); Q_2(2,ny)   = Q_2(nx-4,6); Q_3(2,ny)   = Q_3(nx-4,6); Q_4(2,ny)   = Q_4(nx-4,6)
      Q_1(3,ny-2) = Q_1(nx-3,4); Q_2(3,ny-2) = Q_2(nx-3,4); Q_3(3,ny-2) = Q_3(nx-3,4); Q_4(3,ny-2) = Q_4(nx-3,4)
      Q_1(3,ny-1) = Q_1(nx-3,5); Q_2(3,ny-1) = Q_2(nx-3,5); Q_3(3,ny-1) = Q_3(nx-3,5); Q_4(3,ny-1) = Q_4(nx-3,5)
      Q_1(3,ny)   = Q_1(nx-3,6); Q_2(3,ny)   = Q_2(nx-3,6); Q_3(3,ny)   = Q_3(nx-3,6); Q_4(3,ny)   = Q_4(nx-3,6)
      Q_1(nx-2,ny-2) = Q_1(4,4); Q_2(nx-2,ny-2) = Q_2(4,4); Q_3(nx-2,ny-2) = Q_3(4,4); Q_4(nx-2,ny-2) = Q_4(4,4)
      Q_1(nx-2,ny-1) = Q_1(4,5); Q_2(nx-2,ny-1) = Q_2(4,5); Q_3(nx-2,ny-1) = Q_3(4,5); Q_4(nx-2,ny-1) = Q_4(4,5)
      Q_1(nx-2,ny)   = Q_1(4,6); Q_2(nx-2,ny)   = Q_2(4,6); Q_3(nx-2,ny)   = Q_3(4,6); Q_4(nx-2,ny)   = Q_4(4,6)
      Q_1(nx-1,ny-2) = Q_1(5,4); Q_2(nx-1,ny-2) = Q_2(5,4); Q_3(nx-1,ny-2) = Q_3(5,4); Q_4(nx-1,ny-2) = Q_4(5,4)
      Q_1(nx-1,ny-1) = Q_1(5,5); Q_2(nx-1,ny-1) = Q_2(5,5); Q_3(nx-1,ny-1) = Q_3(5,5); Q_4(nx-1,ny-1) = Q_4(5,5)
      Q_1(nx-1,ny)   = Q_1(5,6); Q_2(nx-1,ny)   = Q_2(5,6); Q_3(nx-1,ny)   = Q_3(5,6); Q_4(nx-1,ny)   = Q_4(5,6)
      Q_1(nx,ny-2)   = Q_1(6,4); Q_2(nx,ny-2)   = Q_2(6,4); Q_3(nx,ny-2)   = Q_3(6,4); Q_4(nx,ny-2)   = Q_4(6,4)
      Q_1(nx,ny-1)   = Q_1(6,5); Q_2(nx,ny-1)   = Q_2(6,5); Q_3(nx,ny-1)   = Q_3(6,5); Q_4(nx,ny-1)   = Q_4(6,5)
      Q_1(nx,ny)     = Q_1(6,6); Q_2(nx,ny)     = Q_2(6,6); Q_3(nx,ny)     = Q_3(6,6); Q_4(nx,ny)     = Q_4(6,6)
    enddo
  end subroutine set_bc_cyclic6
end module set_bc_common
