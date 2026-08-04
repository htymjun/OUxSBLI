program extract
  implicit none
  integer j, nx, nx1, nx2, ny, filesize, ios
  real(8) :: rho, u, v, p, tmp(5)
  real(8) :: gamma = 1.4d0
  real(8), allocatable :: x(:), y(:), Q(:,:,:), Q1d(:,:)
  character(len=40) filename

  write(filename, "(a)") "./recal/x.dat"
  open(10, file=filename, action="read", form="unformatted", access="stream", status="old", iostat=ios) 
  if (ios /= 0) then
    print *, "Error opening file x.dat."
  endif
  inquire(unit=10, size=filesize)
  nx = filesize / 8
  allocate(x(nx))
  read(10) x
  close(10)

  write(filename, "(a)") "./recal/y.dat"
  open(10, file=filename, action="read", form="unformatted", access="stream", status="old", iostat=ios) 
  if (ios /= 0) then
    print *, "Error opening file y.dat."
  endif
  inquire(unit=10, size=filesize)
  ny = filesize / 8
  allocate(y(ny))
  read(10) y
  close(10)

  print *, "Grid size: nx=", nx, ", ny=", ny
  allocate(Q(nx,4,ny))
  
  write(filename, "(a)") "./recal/Q00001.dat"
  open(10, file=filename, action="read", form="unformatted", access="stream", status="old", iostat=ios) 
  if (ios /= 0) then
    print *, "Error opening file Q.dat."
  endif
  read(10) Q
  close(10)

  ! extract surface
  print *, "Choose extract from:"
  read(*,*) nx1
  print *, "to:"
  read(*,*) nx2

  allocate(Q1d(5,ny))
  do j = 1, ny
    rho = sum(Q(nx1:nx2,1,j)) / dble(nx2 - nx1 + 1)
    u   = sum(Q(nx1:nx2,2,j) / Q(nx1:nx2,1,j)) / dble(nx2 - nx1 + 1)
    v   = sum(Q(nx1:nx2,3,j) / Q(nx1:nx2,1,j)) / dble(nx2 - nx1 + 1)
    p   = sum((gamma - 1.d0) * (Q(nx1:nx2,4,j) &
          - 0.5d0 * (Q(nx1:nx2,2,j)**2 + Q(nx1:nx2,3,j)**2) / Q(nx1:nx2,1,j))) / dble(nx2 - nx1 + 1)
    Q1d(1,j) = y(j)
    Q1d(2,j) = rho
    Q1d(3,j) = u
    Q1d(4,j) = v
    Q1d(5,j) = p
  enddo
  ! bc on top
  Q1d(2:5,ny) = Q1d(2:5,ny-1)

  write(filename, "(a, i5.5, a, i5.5, a)") "./Qx", int(nx1), "_", int(nx2), ".d"
  open(10, file=filename, action="write", form="formatted", status="replace") 
  write(10, "(a)") "y rho u v p"
  do j = 1, ny
    tmp = [Q1d(1,j), Q1d(2,j), Q1d(3,j), Q1d(4,j), Q1d(5,j)]
    write(10, *) tmp
  enddo
  close(10)

  write(filename, "(a, i5.5, a, i5.5, a)") "./Qx", int(nx1), "_", int(nx2), ".dat"
  open(10, file=filename, action="write", form="unformatted", access="stream", status="replace") 
  write(10) Q1d
  close(10)

  deallocate(x, y, Q, Q1d)
end program extract
