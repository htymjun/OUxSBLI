#!/usr/bin/env python3
"""CPU check of actual case grid code; CUDA declarations only are stripped."""
from pathlib import Path
import re
import subprocess
import tempfile
case = Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='partition-grid-') as tmp:
    work=Path(tmp)
    (work/'globals.f90').write_text((case/'mod_globals.f90').read_text())
    (work/'grid.f90').write_text((case/'grid_partitioned.f90').read_text())
    s=(case/'set.f90').read_text()
    s=re.sub(r'^\s*use (set_bc_common|set_coordinate)\s*$', '', s, flags=re.M)
    s=re.sub(r',\s*device\b', '', s, flags=re.I)
    s=re.sub(r'attributes\(device\)\s*', '', s, flags=re.I)
    s=s.replace('    call set_bc_mut_common(nx, ny, nz, mut, qc2)', '')
    (work/'set.f90').write_text(s)
    (work/'stub.f90').write_text('''module cudafor
 type dim3
 integer :: x,y,z
 end type
end module
module mod_constant
 integer,parameter :: id_accuracy=0
end module
''')
    (work/'check.f90').write_text('''program check
 use mod_globals
 use set
 use grid_partitioned
 implicit none
 real(8)::x(nx),y(ny),z(nz),dx(nx-1),dy(ny-1),dz(nz-1), xx(nx),yy(ny)
 real(8)::a(nx_main),b(ny_main), h, out(8)
 integer::lo,hi,j
 character(20)::mode
 call get_command_argument(1,mode)
 if(trim(mode)=='invalid') then
 call partitioned_axis(8,1.d0,0.9d0,3,.false.,out)
 stop 9
 endif
 call partitioned_axis(ny_main,Ly_main,Ly_uniform,ny_uniform,.true.,b)
 lo=(ny_main-ny_uniform)/2+1; hi=lo+ny_uniform-1
 h=Ly_uniform/dble(ny_uniform-1)
 call req(maxval(abs(b(lo+1:hi)-b(lo:hi-1)-h))<1.d-13,'y uniform width')
 call req(abs(b(hi+1)-b(hi)-h)<1.d-13,'y junction width')
 call req(maxval(abs(b+b(ny_main:1:-1)-Ly_main))<1.d-12,'y reflection symmetry')
 do j=hi+1,ny_main-1
 call req(b(j+1)-b(j)>=b(j)-b(j-1)-1.d-13,'y outward stretching')
 enddo
 call set_grid(0,nx,ny,nz,Lx,Ly,Lz,x,y,z,dx,dy,dz)
 call req(minval(dx)>0.and.minval(dy)>0.and.minval(dz)>0,'positive full grid widths')
 h=Lx_main/dble(nx_main-1)
 call req(maxval(abs(xi(3:nx_main+1)-xi(2:nx_main)-h))<1.d-13,'entire x main uniform')
 call req(maxval(abs(dx(2:nx_main)-h))<1.d-13,'x main center spacing uniform')
 call req(abs(xi(nx_main+1)-xi(2)-Lx_main)<1.d-12,'x main length')
 call req(abs(xi(nx_main+2)-xi(nx_main+1)-h)<1.d-13,'x buffer junction')
 call req(abs(xi(nx-1)-xi(nx_main+1)-Lx_buf)<1.d-10,'x buffer length')
 call req(abs(yj(ny_buf+2)-yj(2)-Ly_buf)<1.d-10,'lower y buffer length')
 call req(abs(yj(ny-1)-yj(ny_buf+ny_main+1)-Ly_buf)<1.d-10,'upper y buffer length')
 xx=x; yy=y
 call set_grid_main_buffer(nx,ny,nz,Lx_main,0.8d0*Lx_buf,Ly_main,1.2d0*Ly_buf,Lz, &
 nx_main,nx_buf,ny_main,ny_buf,x,y,z,dx,dy,dz)
 call req(maxval(abs(x(2:nx_main)-xx(2:nx_main)))<1.d-12,'x main independent of buffer')
 call req(maxval(abs(y(ny_buf+2:ny_buf+ny_main)-yy(ny_buf+2:ny_buf+ny_main)-0.2d0*Ly_buf)) &
 <1.d-12,'y main widths independent of buffer')
 ! Fully uniform limit is also supported.
 call partitioned_axis(8,1.d0,1.d0,8,.false.,out)
 call req(maxval(abs(out(2:8)-out(1:7)-1.d0/7))<1.d-13,'uniform limit')
 print *, 'PASS: main partitions, symmetry, smooth widths, buffers, independence, uniform limit'
 contains
 subroutine req(ok,label)
 logical,intent(in)::ok
 character(*),intent(in)::label
 if(.not.ok) then
 print *, 'FAIL: ',label
 stop 1
 endif
 end subroutine
end program
''')
    subprocess.run(['gfortran','-O0','-fcheck=all','-fno-range-check','-fwrapv','-ffree-line-length-none',
                    'stub.f90','globals.f90','grid.f90','set.f90','check.f90','-o','check'],cwd=work,check=True)
    subprocess.run([str(work/'check')],cwd=work,check=True)
    invalid=subprocess.run([str(work/'check'),'invalid'],cwd=work,capture_output=True,text=True)
    assert invalid.returncode != 0 and 'Stretch section too short' in invalid.stderr, invalid
    print('PASS: incompatible partition rejected')
