!> VTK output for curvilinear grids (airfoil O-grid)
!> Outputs StructuredGrid format for visualization in ParaView
module print_curv
  use mpi
  use mod_globals, only : step_offset, gamma
  implicit none
  public print_vtk_curv, send_recv_for_print_even_curv, send_recv_for_print_odd_curv, make_1d_for_print_curv
contains
  !> Output solution in VTK StructuredGrid format (for curvilinear grids)
  !> Includes physical coordinates (x_phys, y_phys, z)
  subroutine print_vtk_curv(step, nx, ny, nz, myrank, nranks, &
                             x_phys, y_phys, z, rho_flat, p_flat, vel1d, ke0, entropy0)
    integer, intent(in) :: step, nx, ny, nz, myrank, nranks
    real(8), intent(in) :: x_phys(nx,ny), y_phys(nx,ny), z(nz)
    real(4), intent(in) :: rho_flat(nx*ny*nz), p_flat(nx*ny*nz), vel1d(nx*ny*nz*3)
    real(4), intent(inout) :: ke0, entropy0
    character(len=60) :: filename, result_dir
    integer :: i, j, k, l, m, iunit, ierr
    integer(4) :: byte_x, byte_rho, byte_p, byte_v
    character :: lf*1
    character(len=20) :: extent_str
    real(4), allocatable :: x_flat(:), y_flat(:), z_flat(:)
    real(4), allocatable :: ke_ratio(:), entropy_diff(:)
    
    lf = char(10)
    
    ! Create output directory
    result_dir = 'data'
    call execute_command_line('mkdir -p ' // trim(result_dir), wait=.true., exitstat=ierr)
    
    ! Output filename
    write(filename, '(a, i5.5, a)') trim(result_dir) // '/Q', step+step_offset, '.vts'
    
    ! Allocate flat arrays for StructuredGrid
    allocate(x_flat(nx*ny*nz), y_flat(nx*ny*nz), z_flat(nx*ny*nz))
    allocate(ke_ratio(nx*ny*nz), entropy_diff(nx*ny*nz))
    
    ! Fill flat arrays with 3D coordinates
    ! StructuredGrid expects all nodes in flat format with z varying fastest (Fortran column-major)
    l = 1
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          x_flat(l) = real(x_phys(i,j))
          y_flat(l) = real(y_phys(i,j))
          z_flat(l) = real(z(k))
          ! Compute derived quantities
          ke_ratio(l) = 1.0  ! placeholder
          entropy_diff(l) = 0.0  ! placeholder
          l = l + 1
        enddo
      enddo
    enddo
    
    ! Compute byte sizes for offset calculation.
    ! Each block = 4-byte length header + data bytes.
    ! Points is a 3-component (xyz) array: 4 + 4*3*N bytes total.
    byte_x   = 4 + 4 * 3 * (nx*ny*nz)  ! Points (xyz interleaved)
    byte_rho = 4 + 4 * (nx*ny*nz)
    byte_p   = 4 + 4 * (nx*ny*nz)
    byte_v   = 4 + 4 * (3*nx*ny*nz)
    
    ! Open file for writing
    open(unit=10, file=trim(filename), access='stream', form='unformatted', status='replace')
    
    ! Write VTK XML header and StructuredGrid mesh
    write(extent_str, '(i0, a, i0, a, i0, a, i0, a, i0, a, i0)') &
        0, ' ', nx-1, ' ', 0, ' ', ny-1, ' ', 0, ' ', nz-1
    
    write(10) '<?xml version="1.0"?>' // lf
    write(10) '<VTKFile type="StructuredGrid" version="1.0" byte_order="LittleEndian">' // lf
    write(10) '  <StructuredGrid WholeExtent="' // trim(extent_str) // '">' // lf
    write(10) '    <Piece Extent="' // trim(extent_str) // '">' // lf
    write(10) '      <Points>' // lf
    write(10) '        <DataArray type="Float32" NumberOfComponents="3" format="appended" offset="0"/>' // lf
    write(10) '      </Points>' // lf
    write(10) '      <PointData>' // lf
    write(10) '        <DataArray type="Float32" Name="rho" format="appended" offset="' &
        // trim(int_to_str(byte_x)) // '"/>' // lf
    write(10) '        <DataArray type="Float32" Name="p" format="appended" offset="' &
        // trim(int_to_str(byte_x + byte_rho)) // '"/>' // lf
    write(10) '        <DataArray type="Float32" Name="velocity" NumberOfComponents="3" format="appended" offset="' &
        // trim(int_to_str(byte_x + byte_rho + byte_p)) // '"/>' // lf
    write(10) '      </PointData>' // lf
    write(10) '    </Piece>' // lf
    write(10) '  </StructuredGrid>' // lf
    write(10) '  <AppendedData encoding="raw">' // lf
    write(10) '  _'
    
    ! Write coordinate arrays: one 4-byte length header then all xyz data
    write(10) int(4*3*nx*ny*nz, 4)
    do l = 1, nx*ny*nz
      write(10) x_flat(l), y_flat(l), z_flat(l)
    enddo

    ! Write scalar/vector fields: one 4-byte length header then all data
    write(10) int(4*nx*ny*nz, 4),   rho_flat
    write(10) int(4*nx*ny*nz, 4),   p_flat
    write(10) int(4*3*nx*ny*nz, 4), vel1d
    
    write(10) lf // '  </AppendedData>' // lf
    write(10) '</VTKFile>' // lf
    close(10)
    
    ! Clean up
    deallocate(x_flat, y_flat, z_flat, ke_ratio, entropy_diff)
  end subroutine print_vtk_curv


  !> GPU rank sends solution data; print rank receives and writes VTK.
  subroutine send_recv_for_print_even_curv(myrank, nranks, step, nx, ny, nz, &
      x_phys, y_phys, z, Jacobian_cpu, QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, Q, ke0, entropy0)
    integer, intent(in)         :: myrank, nranks, step, nx, ny, nz
    real(8), intent(in)         :: x_phys(nx,ny), y_phys(nx,ny), z(nz), Jacobian_cpu(nx,ny)
    real(8), intent(in), device :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    real(4), intent(inout)      :: ke0, entropy0
    integer ireq3(3), istat3(MPI_STATUS_SIZE,3), ierr
    real(4) rho_flat(nx*ny*nz), p_flat(nx*ny*nz), vel_flat(nx*ny*nz*3)
    Q(:,:,:,1) = QJ_1
    Q(:,:,:,2) = QJ_2
    Q(:,:,:,3) = QJ_3
    Q(:,:,:,4) = QJ_4
    Q(:,:,:,5) = QJ_5
    call make_1d_for_print_curv(nx, ny, nz, Jacobian_cpu, Q, rho_flat, p_flat, vel_flat)
    call MPI_ISEND(rho_flat, nx*ny*nz,   MPI_REAL4, myrank+1, 3*(myrank+1)-2, MPI_COMM_WORLD, ireq3(1), ierr)
    call MPI_ISEND(p_flat,   nx*ny*nz,   MPI_REAL4, myrank+1, 3*(myrank+1)-1, MPI_COMM_WORLD, ireq3(2), ierr)
    call MPI_ISEND(vel_flat,   nx*ny*nz*3, MPI_REAL4, myrank+1, 3*(myrank+1),   MPI_COMM_WORLD, ireq3(3), ierr)
    call MPI_WAITALL(3, ireq3, istat3, ierr)
  end subroutine send_recv_for_print_even_curv


  subroutine send_recv_for_print_odd_curv(myrank, nranks, step, nx, ny, nz, &
      x_phys, y_phys, z, Jacobian_cpu, Q, ke0, entropy0)
    integer, intent(in)    :: myrank, nranks, step, nx, ny, nz
    real(8), intent(in)    :: x_phys(nx,ny), y_phys(nx,ny), z(nz), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    real(4), intent(inout) :: ke0, entropy0
    integer ireq3(3), istat3(MPI_STATUS_SIZE,3), ierr
    real(4) rho_flat(nx*ny*nz), p_flat(nx*ny*nz), vel_flat(nx*ny*nz*3)
    call MPI_IRECV(rho_flat, nx*ny*nz,   MPI_REAL4, myrank-1, 3*myrank-2, MPI_COMM_WORLD, ireq3(1), ierr)
    call MPI_IRECV(p_flat,   nx*ny*nz,   MPI_REAL4, myrank-1, 3*myrank-1, MPI_COMM_WORLD, ireq3(2), ierr)
    call MPI_IRECV(vel_flat,   nx*ny*nz*3, MPI_REAL4, myrank-1, 3*myrank,   MPI_COMM_WORLD, ireq3(3), ierr)
    call MPI_WAITALL(3, ireq3, istat3, ierr)
    call print_vtk_curv(step, nx, ny, nz, myrank, nranks, x_phys, y_phys, z, rho_flat, p_flat, vel_flat, ke0, entropy0)
  end subroutine send_recv_for_print_odd_curv


  !> Helper function to convert integer to string for offset calculation
  function int_to_str(i) result(s)
    integer(4), intent(in) :: i
    character(len=20) :: s
    write(s, '(i0)') i
    s = adjustl(s)
  end function int_to_str


  !> Convert curvilinear conserved variable QJ to real(4) 1D arrays for VTK output.
  !> QJ(nx,5,ny,nz) = [rho/J, rho*u/J, rho*v/J, rho*w/J, rho*E/J]; Jacobian = 1/(J_2D*dz).
  subroutine make_1d_for_print_curv(nx, ny, nz, Jacobian, QJ, rho_flat, p_flat, vel_flat)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: Jacobian(nx,ny), QJ(nx,ny,nz,5)
    real(4), intent(out), dimension(nx*ny*nz)   :: rho_flat, p_flat
    real(4), intent(out), dimension(nx*ny*nz*3) :: vel_flat
    real(8) rho, u, v, w, p
    integer i, j, k, l, m
    l = 1
    m = 1
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          rho      = Jacobian(i,j) * QJ(i,j,k,1)
          u        = QJ(i,j,k,2) / QJ(i,j,k,1)
          v        = QJ(i,j,k,3) / QJ(i,j,k,1)
          w        = QJ(i,j,k,4) / QJ(i,j,k,1)
          p        = (gamma - 1.d0) * (Jacobian(i,j) * QJ(i,j,k,5) - 0.5d0 * rho * (u**2 + v**2 + w**2))
          rho_flat(l) = real(rho, 4)
          p_flat(l)   = real(p,   4)
          vel_flat(m)   = real(u,   4)
          vel_flat(m+1) = real(v,   4)
          vel_flat(m+2) = real(w,   4)
          l = l + 1
          m = m + 3
    enddo;enddo;enddo
  end subroutine make_1d_for_print_curv
end module print_curv
