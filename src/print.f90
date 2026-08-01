#define io 4
module print
  use mpi
  use mod_globals, only : nt, np, dt, step_offset, gamma, R, dimension
  use mod_constant, only : id_accuracy
  implicit none
  integer, parameter :: IO_UNIT_VTK = 10
  interface
    subroutine print_entropy(step, nx, ny, nz, rho_flat, p_flat, entropy0, myrank)
      integer, intent(in)           :: step, nx, ny, nz
      real(io), intent(in)          :: rho_flat(nx*ny*nz), p_flat(nx*ny*nz)
      real(4), intent(inout)        :: entropy0
      integer, intent(in), optional :: myrank
    end subroutine print_entropy
  end interface
 
  interface
    subroutine print_KE(step, nx, ny, nz, rho_flat, vel_flat, ke0, myrank)
      integer, intent(in)           :: step, nx, ny, nz
      real(io), intent(in)          :: rho_flat(nx*ny*nz), vel_flat(nx*ny*nz*3)
      real(4), intent(inout)        :: ke0
      integer, intent(in), optional :: myrank
    end subroutine print_KE
  end interface

  interface make_1d_for_print
    module procedure make_1d_for_print2, make_1d_for_print3
  end interface
  
  interface send_recv_for_print_even
    module procedure send_recv_for_print_even2, send_recv_for_print_even3
  end interface

  interface send_recv_for_print_odd
    module procedure send_recv_for_print_odd2, send_recv_for_print_odd3
  end interface

  interface print_vtk
    module procedure print_vtk_2D, print_vtk_3D
  end interface
contains
  subroutine print_entropy(step, nx, ny, nz, rho_flat, p_flat, entropy0, myrank)
    integer, intent(in)           :: step, nx, ny, nz
    real(io), intent(in)          :: rho_flat(nx*ny*nz), p_flat(nx*ny*nz)
    real(4), intent(inout)        :: entropy0
    integer, intent(in), optional :: myrank
    real(io) entropy, t
    character(len=40) filename
    integer i, j, k, l, accuracy, offset
    if (kind(id_accuracy) == 8) then
      accuracy = 6
      offset   = accuracy / 2
    elseif (kind(id_accuracy) == 4) then
      accuracy = 4
      offset   = accuracy / 2
    else
      accuracy = 2
      offset   = accuracy / 2
    endif
    entropy = 0.e0
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          l = i + (j-1) * nx + (k-1) * nx * ny
          entropy = entropy + rho_flat(l) * log(p_flat(l) * (rho_flat(l)**(-gamma)))
    enddo;enddo;enddo
    entropy = entropy / dble((nx-accuracy) * (ny-accuracy) * (nz-accuracy))
    if (step == 0) then
      entropy0 = entropy
    endif
    t = nt * step * dt
    if (present(myrank)) then
      write(filename, "(a, i0, a)") "data/",int(myrank),"/entropy.d"
      open(IO_UNIT_VTK,file=filename, position="append")
    else
      open(IO_UNIT_VTK,file="data/entropy.d", position="append")
    endif
    write(IO_UNIT_VTK,"(2e12.4)") t, (entropy0 - entropy) / entropy0
    close(IO_UNIT_VTK)
  end subroutine print_entropy


  subroutine print_KE(step, nx, ny, nz, rho_flat, vel_flat, ke0, myrank)
    integer, intent(in)           :: step, nx, ny, nz
    real(io), intent(in)          :: rho_flat(nx*ny*nz), vel_flat(nx*ny*nz*3)
    real(4), intent(inout)        :: ke0
    integer, intent(in), optional :: myrank
    real(io) ke, t
    character(len=40) filename
    integer i, j, k, l, m, accuracy, offset
    if (kind(id_accuracy) == 8) then
      accuracy = 6
      offset   = accuracy / 2
    elseif (kind(id_accuracy) == 4) then
      accuracy = 4
      offset   = accuracy / 2
    else
      accuracy = 2
      offset   = accuracy / 2
    endif
    ke = 0.e0
    do k = 1+offset, nz-offset
      do j = 1+offset, ny-offset
        do i = 1+offset, nx-offset
          l  = i + (j-1) * nx + (k-1) * nx * ny
          m  = 1 + (i-1) * 3 + (j-1) * 3 * nx + (k-1) * 3 * nx * ny 
          ke = ke + 0.5d0 * rho_flat(l) * (vel_flat(m)**2 + vel_flat(m+1)**2 + vel_flat(m+2)**2)
    enddo;enddo;enddo
    ke = ke / dble((nx-accuracy) * (ny-accuracy) * (nz-accuracy))
    if (step == 0) then
      ke0 = ke
    endif
    t = nt * step * dt
    if (present(myrank)) then
      write(filename, "(a, i0, a)") "data/",int(myrank),"/kinetic_energy.d"
      open(IO_UNIT_VTK,file=filename, position="append")
    else
      open(IO_UNIT_VTK,file="data/kinetic_energy.d", position="append")
    endif
    write(IO_UNIT_VTK,"(3e12.4)") t, ke, ke / ke0
    close(IO_UNIT_VTK)
  end subroutine print_KE


  subroutine print_xml(ni, nj, nk, dimension, x, y, z, rho_flat, p_flat, vel_flat)
    integer, intent(in)  :: ni, nj, nk, dimension
    real(io), intent(in) :: x(ni), y(nj), z(nk)
    real(io), intent(in) :: rho_flat(ni*nj*nk), p_flat(ni*nj*nk)
    real(io), intent(in) :: vel_flat(dimension*ni*nj*nk)
    integer(4) byte_x, byte_y, byte_z, byte_rho, byte_p, byte_v
    character :: lf*1, str1*4, str2*4, str3*4, str4*1
    character :: offset1*12, offset2*12, offset3*12, offset4*12, offset5*12
    lf = char(10)
    write(str1(1:4),'(i4)') ni-1
    write(str2(1:4),'(i4)') nj-1
    write(str3(1:4),'(i4)') nk-1
    write(str4(1:1),'(i1)') dimension
    byte_x   = 4 + io * ni
    byte_y   = 4 + io * nj
    byte_z   = 4 + io * nk
    byte_rho = 4 + io * (ni * nj * nk)
    byte_p   = byte_rho
    byte_v   = 4 + io * (dimension * ni * nj * nk)
    write(offset1(1:12),'(i12)') int(byte_x, kind=8)
    write(offset2(1:12),'(i12)') int(byte_x, kind=8) + int(byte_y, kind=8)
    write(offset3(1:12),'(i12)') int(byte_x, kind=8) + int(byte_y, kind=8) + int(byte_z, kind=8)
    write(offset4(1:12),'(i12)') int(byte_x, kind=8) + int(byte_y, kind=8) + int(byte_z, kind=8) + &
                                 int(byte_rho, kind=8)
    write(offset5(1:12),'(i12)') int(byte_x, kind=8) + int(byte_y, kind=8) + int(byte_z, kind=8) + &
                                 int(byte_rho, kind=8) + int(byte_p, kind=8)
    write(IO_UNIT_VTK) '<?xml version="1.0"?>'//lf
    write(IO_UNIT_VTK) '<VTKFile type="RectilinearGrid" version="1.0" byte_order="LittleEndian">'//lf
    write(IO_UNIT_VTK) '  <RectilinearGrid WholeExtent="0 '//str1//' 0 '//str2//' 0 '//str3//'">'//lf
    write(IO_UNIT_VTK) '    <Piece Extent="0 '//str1//' 0 '//str2//' 0 '//str3//'">'//lf
    write(IO_UNIT_VTK) '      <Coordinates>'//lf
    if (io == 4) then
      write(IO_UNIT_VTK) '      <DataArray type="Float32" format="appended" offset="0"/>'//lf
      write(IO_UNIT_VTK) '      <DataArray type="Float32" format="appended" offset="'//offset1//'"/>'//lf
      write(IO_UNIT_VTK) '      <DataArray type="Float32" format="appended" offset="'//offset2//'"/>'//lf
      write(IO_UNIT_VTK) '    </Coordinates>'//lf
      write(IO_UNIT_VTK) '    <PointData>'//lf
      write(IO_UNIT_VTK) '      <DataArray type="Float32" format="appended" offset="'//offset3//'"&
      & Name="rho" NumberOfComponents="1"/>'//lf
      write(IO_UNIT_VTK) '      <DataArray type="Float32" format="appended" offset="'//offset4//'"&
      & Name="p" NumberOfComponents="1"/>'//lf
      write(IO_UNIT_VTK) '      <DataArray type="Float32" format="appended" offset="'//offset5//'"&
      & Name="velocity" NumberOfComponents="'//str4//'"/>'//lf
    else
      write(IO_UNIT_VTK) '      <DataArray type="Float64" format="appended" offset="0"/>'//lf
      write(IO_UNIT_VTK) '      <DataArray type="Float64" format="appended" offset="'//offset1//'"/>'//lf
      write(IO_UNIT_VTK) '      <DataArray type="Float64" format="appended" offset="'//offset2//'"/>'//lf
      write(IO_UNIT_VTK) '    </Coordinates>'//lf
      write(IO_UNIT_VTK) '    <PointData>'//lf
      write(IO_UNIT_VTK) '      <DataArray type="Float64" format="appended" offset="'//offset3//'"&
      & Name="rho" NumberOfComponents="1"/>'//lf
      write(IO_UNIT_VTK) '      <DataArray type="Float64" format="appended" offset="'//offset4//'"&
      & Name="p" NumberOfComponents="1"/>'//lf
      write(IO_UNIT_VTK) '      <DataArray type="Float64" format="appended" offset="'//offset5//'"&
      & Name="velocity" NumberOfComponents="'//str4//'"/>'//lf
    endif
    write(IO_UNIT_VTK) '      </PointData>'//lf
    write(IO_UNIT_VTK) '    </Piece>'//lf
    write(IO_UNIT_VTK) '  </RectilinearGrid>'//lf
    write(IO_UNIT_VTK) '  <AppendedData encoding="raw">'//lf
    write(IO_UNIT_VTK) '  _', byte_x, x, byte_y, y, byte_z, z, byte_rho, rho_flat, byte_p, p_flat, byte_v, vel_flat, lf
    write(IO_UNIT_VTK) '  </AppendedData>'//lf
    write(IO_UNIT_VTK) '</VTKFile>'//lf
    close(IO_UNIT_VTK)
  end subroutine print_xml


  subroutine make_1d_for_print2(nx, ny, Jacobian, QJ, rho_flat, p_flat, vel_flat)
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: Jacobian(nx,ny), QJ(nx,4,ny) ! Q / Jacobian
    real(io), intent(out), dimension(nx*ny)   :: rho_flat, p_flat
    real(io), intent(out), dimension(nx*ny*3) :: vel_flat
    real(8) rho, u, v, p
    integer i, j, k, l, m
    l = 1
    m = 1
    do j = 1, ny
      do i = 1, nx
        rho      = Jacobian(i,j) * QJ(i,1,j)
        u        = QJ(i,2,j) / QJ(i,1,j)
        v        = QJ(i,3,j) / QJ(i,1,j)
        p        = (gamma - 1.d0) * (Jacobian(i,j) * QJ(i,4,j) - 0.5d0 * rho * (u**2 + v**2))
        rho_flat(l) = real(rho, kind=io)
        p_flat(l)   = real(p,   kind=io)
        vel_flat(m)   = real(u,   kind=io)
        vel_flat(m+1) = real(v,   kind=io)
        vel_flat(m+2) = 0.e0
        l = l + 1
        m = m + 3
    enddo;enddo
  end subroutine make_1d_for_print2


  subroutine make_1d_for_print3(nx, ny, nz, Jacobian, QJ, rho_flat, p_flat, vel_flat)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: Jacobian(nx,ny), QJ(nx,ny,nz,5) ! Q / Jacobian
    real(io), intent(out), dimension(nx*ny*nz)   :: rho_flat, p_flat
    real(io), intent(out), dimension(nx*ny*nz*3) :: vel_flat
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
          rho_flat(l) = real(rho, kind=io)
          p_flat(l)   = real(p,   kind=io)
          vel_flat(m)   = real(u,   kind=io)
          vel_flat(m+1) = real(v,   kind=io)
          vel_flat(m+2) = real(w,   kind=io)
          l = l + 1
          m = m + 3
    enddo;enddo;enddo
  end subroutine make_1d_for_print3


  subroutine send_recv_for_print_even2(myrank, nranks, step, nx, ny, x, y, Jacobian_cpu, QJ, Q, ke0, entropy0)
    integer, intent(in)         :: myrank, nranks, step, nx, ny
    real(8), intent(in)         :: x(nx), y(ny), Jacobian_cpu(nx,ny)
    real(8), intent(in), device :: QJ(nx,4,ny)
    real(8), intent(inout)      :: Q(nx,4,ny)
    real(4), intent(inout)      :: ke0, entropy0
    integer ireq3(3), istat3(MPI_STATUS_SIZE,3), ierr
    real(io) rho_flat(nx*ny), p_flat(nx*ny), vel_flat(nx*ny*3)
    Q = QJ
    call make_1d_for_print(nx, ny, Jacobian_cpu, Q, rho_flat, p_flat, vel_flat)
    if (io == 4) then
      call MPI_ISEND(rho_flat, nx*ny,   MPI_REAL4, myrank+1, 3*(myrank+1)-2, MPI_COMM_WORLD, ireq3(1), ierr)
      call MPI_ISEND(p_flat,   nx*ny,   MPI_REAL4, myrank+1, 3*(myrank+1)-1, MPI_COMM_WORLD, ireq3(2), ierr)
      call MPI_ISEND(vel_flat,   nx*ny*3, MPI_REAL4, myrank+1, 3*(myrank+1),   MPI_COMM_WORLD, ireq3(3), ierr)
    else
      call MPI_ISEND(rho_flat, nx*ny,   MPI_REAL8, myrank+1, 3*(myrank+1)-2, MPI_COMM_WORLD, ireq3(1), ierr)
      call MPI_ISEND(p_flat,   nx*ny,   MPI_REAL8, myrank+1, 3*(myrank+1)-1, MPI_COMM_WORLD, ireq3(2), ierr)
      call MPI_ISEND(vel_flat,   nx*ny*3, MPI_REAL8, myrank+1, 3*(myrank+1),   MPI_COMM_WORLD, ireq3(3), ierr)
    endif
    call MPI_WAITALL(3, ireq3, istat3, ierr)
  end subroutine send_recv_for_print_even2


  subroutine send_recv_for_print_even3(myrank, nranks, step, nx, ny, nz, x, y, z, Jacobian_cpu, &
                                       QJ_1, QJ_2, QJ_3, QJ_4, QJ_5, Q, ke0, entropy0)
    integer, intent(in)         :: myrank, nranks, step, nx, ny, nz
    real(8), intent(in)         :: x(nx), y(ny), z(nz), Jacobian_cpu(nx,ny)
    real(8), intent(in), device :: QJ_1(nx,ny,nz), QJ_2(nx,ny,nz), QJ_3(nx,ny,nz), QJ_4(nx,ny,nz), QJ_5(nx,ny,nz)
    real(8), intent(inout)      :: Q(nx,ny,nz,5)
    real(4), intent(inout)      :: ke0, entropy0
    integer ireq3(3), istat3(MPI_STATUS_SIZE,3), ierr
    real(io) rho_flat(nx*ny*nz), p_flat(nx*ny*nz), vel_flat(nx*ny*nz*3)
    Q(:,:,:,1) = QJ_1
    Q(:,:,:,2) = QJ_2
    Q(:,:,:,3) = QJ_3
    Q(:,:,:,4) = QJ_4
    Q(:,:,:,5) = QJ_5
    call make_1d_for_print(nx, ny, nz, Jacobian_cpu, Q, rho_flat, p_flat, vel_flat)
    if (io == 4) then
      call MPI_ISEND(rho_flat, nx*ny*nz,   MPI_REAL4, myrank+1, 3*(myrank+1)-2, MPI_COMM_WORLD, ireq3(1), ierr)
      call MPI_ISEND(p_flat,   nx*ny*nz,   MPI_REAL4, myrank+1, 3*(myrank+1)-1, MPI_COMM_WORLD, ireq3(2), ierr)
      call MPI_ISEND(vel_flat,   nx*ny*nz*3, MPI_REAL4, myrank+1, 3*(myrank+1),   MPI_COMM_WORLD, ireq3(3), ierr)
    else
      call MPI_ISEND(rho_flat, nx*ny*nz,   MPI_REAL8, myrank+1, 3*(myrank+1)-2, MPI_COMM_WORLD, ireq3(1), ierr)
      call MPI_ISEND(p_flat,   nx*ny*nz,   MPI_REAL8, myrank+1, 3*(myrank+1)-1, MPI_COMM_WORLD, ireq3(2), ierr)
      call MPI_ISEND(vel_flat,   nx*ny*nz*3, MPI_REAL8, myrank+1, 3*(myrank+1),   MPI_COMM_WORLD, ireq3(3), ierr)
    endif
    call MPI_WAITALL(3, ireq3, istat3, ierr)
  end subroutine send_recv_for_print_even3


  subroutine send_recv_for_print_odd2(myrank, nranks, step, nx, ny, x, y, Jacobian_cpu, Q, ke0, entropy0)
    integer, intent(in)    :: myrank, nranks, step, nx, ny
    real(8), intent(in)    :: x(nx), y(ny), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(nx,4,ny)
    real(4), intent(inout) :: ke0, entropy0
    integer ireq3(3), istat3(MPI_STATUS_SIZE,3), ierr
    real(io) rho_flat(nx*ny), p_flat(nx*ny), vel_flat(nx*ny*3)
    if (io == 4) then
      call MPI_IRECV(rho_flat, nx*ny,   MPI_REAL4, myrank-1, 3*myrank-2, MPI_COMM_WORLD, ireq3(1), ierr)
      call MPI_IRECV(p_flat,   nx*ny,   MPI_REAL4, myrank-1, 3*myrank-1, MPI_COMM_WORLD, ireq3(2), ierr)
      call MPI_IRECV(vel_flat,   nx*ny*3, MPI_REAL4, myrank-1, 3*myrank,   MPI_COMM_WORLD, ireq3(3), ierr)
    else
      call MPI_IRECV(rho_flat, nx*ny,   MPI_REAL8, myrank-1, 3*myrank-2, MPI_COMM_WORLD, ireq3(1), ierr)
      call MPI_IRECV(p_flat,   nx*ny,   MPI_REAL8, myrank-1, 3*myrank-1, MPI_COMM_WORLD, ireq3(2), ierr)
      call MPI_IRECV(vel_flat,   nx*ny*3, MPI_REAL8, myrank-1, 3*myrank,   MPI_COMM_WORLD, ireq3(3), ierr)
    endif
    call MPI_WAITALL(3, ireq3, istat3, ierr)
    call print_vtk(step, nx, ny, x, y, rho_flat, p_flat, vel_flat)
  end subroutine send_recv_for_print_odd2


  subroutine send_recv_for_print_odd3(myrank, nranks, step, nx, ny, nz, x, y, z, Jacobian_cpu, Q, ke0, entropy0)
    integer, intent(in)    :: myrank, nranks, step, nx, ny, nz
    real(8), intent(in)    :: x(nx), y(ny), z(nz), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(nx,ny,nz,5)
    real(4), intent(inout) :: ke0, entropy0
    integer ireq3(3), istat3(MPI_STATUS_SIZE,3), ierr
    real(io) rho_flat(nx*ny*nz), p_flat(nx*ny*nz), vel_flat(nx*ny*nz*3)
    if (io == 4) then
      call MPI_IRECV(rho_flat, nx*ny*nz,   MPI_REAL4, myrank-1, 3*myrank-2, MPI_COMM_WORLD, ireq3(1), ierr)
      call MPI_IRECV(p_flat,   nx*ny*nz,   MPI_REAL4, myrank-1, 3*myrank-1, MPI_COMM_WORLD, ireq3(2), ierr)
      call MPI_IRECV(vel_flat,   nx*ny*nz*3, MPI_REAL4, myrank-1, 3*myrank,   MPI_COMM_WORLD, ireq3(3), ierr)
    else
      call MPI_IRECV(rho_flat, nx*ny*nz,   MPI_REAL8, myrank-1, 3*myrank-2, MPI_COMM_WORLD, ireq3(1), ierr)
      call MPI_IRECV(p_flat,   nx*ny*nz,   MPI_REAL8, myrank-1, 3*myrank-1, MPI_COMM_WORLD, ireq3(2), ierr)
      call MPI_IRECV(vel_flat,   nx*ny*nz*3, MPI_REAL8, myrank-1, 3*myrank,   MPI_COMM_WORLD, ireq3(3), ierr)
    endif
    call MPI_WAITALL(3, ireq3, istat3, ierr)
    call print_vtk(step, nx, ny, nz, myrank, nranks, x, y, z, rho_flat, p_flat, vel_flat, ke0, entropy0)
  end subroutine send_recv_for_print_odd3


  subroutine print0(myrank, nx, ny, x, y, Jacobian_cpu, Q, ke0, entropy0)
    integer, intent(in)    :: myrank, nx, ny
    real(8), intent(in)    :: x(nx)               !< x coordinate array (host)
    real(8), intent(in)    :: y(ny)               !< y coordinate array (host)
    real(8), intent(in)    :: Jacobian_cpu(nx,ny) !< Jacobian determinant (host)
    real(8), intent(in)    :: Q(nx,4,ny)          !< conservative variables on host
    real(4), intent(inout) :: ke0
    real(4), intent(inout) :: entropy0
    real(io) rho_flat(nx*ny), p_flat(nx*ny), vel_flat(nx*ny*3)
    integer ierr
    call make_1d_for_print(nx, ny, Jacobian_cpu, Q, rho_flat, p_flat, vel_flat)
    call print_vtk(0, nx, ny, x, y, rho_flat, p_flat, vel_flat)
    call MPI_SEND(ke0,      1, MPI_REAL4, myrank+1, 3*(myrank+1)+1, MPI_COMM_WORLD, ierr)
    call MPI_SEND(entropy0, 1, MPI_REAL4, myrank+1, 3*(myrank+1)+2, MPI_COMM_WORLD, ierr)
  end subroutine print0


  subroutine print_vtk_2D(step, nx, ny, x, y, rho_flat, p_flat, vel_flat)
    integer, intent(in)  :: step, nx, ny
    real(8), intent(in)  :: x(nx), y(ny)
    real(io), intent(in) :: rho_flat(nx*ny), p_flat(nx*ny), vel_flat(nx*ny*3)
    real(8) :: z(1) = 0.d0
    character(len=40) filename
    write(filename, "(a, i5.5,a)") "data/Q",int(step+step_offset),".vtr"
    open(IO_UNIT_VTK,file=filename,status="replace",action="write",form="unformatted",access="stream",convert="Little_ENDIAN")
    call print_xml(nx, ny, 1, 3, real(x, kind=io), real(y, kind=io), real(z, kind=io), rho_flat, p_flat, vel_flat)
    close(IO_UNIT_VTK)
  end subroutine print_vtk_2D


  subroutine print_vtk_3D(step, nx, ny, nz, myrank, nranks, x, y, z, rho_flat, p_flat, vel_flat, ke0, entropy0)
    integer, intent(in)    :: step, nx, ny, nz, myrank, nranks
    real(8), intent(in)    :: x(nx), y(ny), z(nz)
    real(io), intent(in)   :: rho_flat(nx*ny*nz), p_flat(nx*ny*nz), vel_flat(nx*ny*nz*3)
    real(4), intent(inout) :: ke0, entropy0
    character(len=40) filename
    if (nranks >= 4) then
      call print_entropy(step, nx, ny, nz, rho_flat, p_flat, entropy0, myrank)
      call print_KE(step, nx, ny, nz, rho_flat, vel_flat, ke0, myrank)
      write(filename, "(a, i0, a, i5.5, a)") "data/",int(myrank),"/Q",int(step+step_offset),".vtr"
    else
      call print_entropy(step, nx, ny, nz, rho_flat, p_flat, entropy0)
      call print_KE(step, nx, ny, nz, rho_flat, vel_flat, ke0)
      write(filename, "(a, i5.5, a)") "data/Q",int(step+step_offset),".vtr"
    endif
    open(IO_UNIT_VTK,file=filename,status="replace",action="write",form="unformatted",access="stream",convert="Little_ENDIAN")
    call print_xml(nx, ny, nz, 3, real(x, kind=io), real(y, kind=io), real(z, kind=io), rho_flat, p_flat, vel_flat)
    close(IO_UNIT_VTK)
  end subroutine print_vtk_3D
end module print
