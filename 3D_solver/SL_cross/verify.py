#!/usr/bin/env python3
"""Run the case's Fortran initialization/BC loops on the CPU at reduced size.

CUDA attributes/directives and the x mesh size are adapted in temporary files.
Checks case logic and the shared residual routine, but not GPU execution.
Requires gfortran; no GPU or third-party Python packages are needed.
"""
from pathlib import Path
import re
import subprocess
import tempfile

case = Path(__file__).resolve().parent
with tempfile.TemporaryDirectory(prefix='sl-cross-check-') as directory:
    work = Path(directory)
    globals_text = (case / 'mod_globals.f90').read_text()
    globals_text = re.sub(r'(:: nx_main =)[^\n]+', r'\1 16', globals_text)
    globals_text = globals_text.replace('inlet_fluctuation_rms = 0.05d0', 'inlet_fluctuation_rms = 0.d0')
    (work / 'mod_globals.f90').write_text(globals_text)
    source = (case / 'set.f90').read_text()
    source = re.sub(r'^\s*use (set_bc_common|set_coordinate)\s*$', '', source, flags=re.M)
    source = re.sub(r',\s*device\b', '', source, flags=re.I)
    source = re.sub(r'attributes\(device\)\s*', '', source, flags=re.I)
    (work / 'set.f90').write_text(source)
    # Exercise the actual shared residual routine with and without z correction.
    steps = (case.parent / 'src/calc_steps.f90').read_text()
    steps = steps[:steps.index('  end subroutine calc_R')]+ '  end subroutine calc_R\nend module calc_steps\n'
    steps = re.sub(r',\s*(device|contiguous)\b', '', steps, flags=re.I)
    steps = re.sub(r'attributes\(device\)\s*', '', steps, flags=re.I)
    (work / 'calc_steps.F90').write_text(steps)
    (work / 'stubs.f90').write_text('''module cudafor
 type dim3
  integer :: x, y, z
 end type
end module
module mod_constant
 integer, parameter :: id_accuracy=0
 real(8), parameter :: one_sixth=1.d0/6.d0
end module
''')
    (work / 'check.f90').write_text('''program check
 use mod_globals
 use set
 use calc_steps, only : calc_R
 implicit none
 real(8) :: x(nx),y(ny),z(nz),dx(nx-1),dy(ny-1),dz(nz-1),jac(nx,ny)
 real(8) :: q(nx,5,ny,nz),before(nx,5,ny,nz), a,b,t,density_value,e, expected, err
 real(8) :: faces(nz+1), ef(5,3,2,nz-2),ff(5,2,3,nz-2),gf(5,2,2,nz-1),res(5)
 real(8) :: hx,hy,hz,dzref, weights(nz), total, exact
 integer :: i,j,k,l,axis
 call set_grid(0,nx,ny,nz,Lx,Ly,Lz,x,y,z,dx,dy,dz)
 call require(minval(dx)>0.and.minval(dy)>0.and.minval(dz)>0, 'positive spacing')
 call require(ny==nz,'matching transverse grid counts')
 call require(maxval(abs(y-z))<1.d-13,'identical y/z coordinates')
 call require(maxval(dz)>2*minval(dz),'nonuniform z mesh')
 call transverse_faces(nz,nz_main,nz_buf,Lz_main,Lz_buf,faces)
 call require(abs(faces(nz_buf+2)-faces(2)-Lz_buf)<1.d-10,'lower z buffer length')
 call require(abs(faces(nz-1)-faces(nz_buf+nz_main+1)-Lz_buf)<1.d-10,'upper z buffer length')
 call require(abs(faces(nz_buf+nz_main+1)-faces(nz_buf+2)-Lz_main)<1.d-12,'z core length')
 call require(maxval(abs(sigma_z_1d-sigma_z_1d(nz:1:-1)))<1.d-14, 'symmetric z sponge')
 call require(maxval(sigma_z_1d(nz_buf+2:nz_buf+nz_main+1))==0.d0, 'no damping in z core')
 call cross_profile(Ly/2+5*delta_bl,Lz/2+5*delta_bl,a,t,density_value,e)
 call cross_profile(Ly/2+5*delta_bl,Lz/2-5*delta_bl,b,t,density_value,e)
 call require(abs(a-u1)<1.d-5.and.abs(b-u2)<1.d-5, 'high/low quadrants')
 call cross_profile(Ly/2-5*delta_bl,Lz/2-5*delta_bl,a,t,density_value,e)
 call cross_profile(Ly/2-5*delta_bl,Lz/2+5*delta_bl,b,t,density_value,e)
 call require(abs(a-u1)<1.d-5.and.abs(b-u2)<1.d-5, 'reversed quadrants')
 call cross_profile(Ly/2,Lz/2,a,t,density_value,e)
 call require(abs(a-(u1+u2)/2)<1.d-12.and.abs(e-1)<1.d-12, 'cross center')
 call set_Jacobian_xy3_stretch(nx,ny,nz,x,y,z,jac)
 call set_init(0,nx,ny,nz,x,y,z,q)
 do k=1,nz
  do j=1,ny
   call require(abs(q(2,2,j,k)/q(2,1,j,k)-u_target(j,k))<1.d-10,'cached velocity')
   call require(abs(q(2,1,j,k)-rho_target(j,k))<1.d-12,'cached density')
   expected=(gamma-1)*(q(2,5,j,k)-0.5d0*q(2,2,j,k)**2/q(2,1,j,k))
   call require(abs(expected-p)<1.d-7,'initial pressure')
   do l=1,5
    q(:,l,j,k)=q(:,l,j,k)/jac(:,j)
   enddo
  enddo
 enddo
 before=q
 call set_bc(0,nx,ny,nz,jac,q)
 err=maxval(abs(q(2:nx-1,:,2:ny-1,2:nz-1)-before(2:nx-1,:,2:ny-1,2:nz-1)))
 call require(err<1.d-10,'sponge preserves unperturbed cross field')
 do k=2,nz-1
  do j=2,ny-1
   call require(maxval(abs(q(1,:,j,k)-before(1,:,j,k)))<1.d-10,'inlet matches initial field')
  enddo
 enddo
 sigma_x_1d=0; sigma_y_1d=0; sigma_z_1d=0
 do k=1,nz
  do j=1,ny
   do i=1,nx
    do l=1,5
     q(i,l,j,k)=(10.d0*l+0.01d0*i+0.02d0*j+0.03d0*z(k))/jac(i,j)
    enddo
   enddo
  enddo
 enddo
 call set_bc(0,nx,ny,nz,jac,q)
 do k=1,nz
  do j=1,ny
   do i=2,nx
    do l=1,5
     expected=10.d0*l+0.01d0*i+0.02d0*max(2,min(ny-1,j))+0.03d0*z(max(2,min(nz-1,k)))
     call require(abs(q(i,l,j,k)*jac(i,j)-expected)<1.d-10,'x linear, yz constant extrapolation incl corners')
    enddo
   enddo
  enddo
 enddo
 ! Manufactured linear face fluxes: each physical divergence must be one.
 hx=0.5d0*(x(4)-x(2)); hy=0.5d0*(y(4)-y(2))
 dzref=0.5d0*(z(nz/2+1)-z(nz/2-1))
 weights=z_residual_scale
 do axis=1,3
  do k=1,nz-2
   hz=0.5d0*(z(k+2)-z(k))
   ef=0; ff=0; gf=0
   if(axis==1) ef(:,2,1,k)=hx
   if(axis==2) ff(:,1,2,k)=hy
   if(axis==3) gf(:,1,1,k+1)=hz
   call calc_R(4,4,nz,1,1,k,dt*hx*hy,dt*hy*hz,dt*hx*hz,ef,ff,gf,res)
   expected=dt*hx*hy*dzref
   call require(maxval(abs(res/expected-1.d0))<1.d-12,'xyz metric manufactured divergence')
  enddo
 enddo
 ! Physical-volume weighted sum must telescope to the outer z fluxes.
 ef=0; ff=0
 do k=1,nz-1
  gf(:,:,:,k)=sin(0.07d0*k)
 enddo
 total=0
 do k=1,nz-2
  hz=0.5d0*(z(k+2)-z(k))
  call calc_R(4,4,nz,1,1,k,dt*hx*hy,dt*hy*hz,dt*hx*hz,ef,ff,gf,res)
  total=total+res(1)/weights(k+1)
 enddo
 exact=dt*hx*hy*(gf(1,1,1,nz-1)-gf(1,1,1,1))
 call require(abs(total-exact)<1.d-12*dt*hx*hy,'volume-weighted conservation')
 ! Constant flux gives zero; uniform z reduces to the original residual.
 ef=1; ff=1; gf=1
 call calc_R(4,4,nz,1,1,1,dt,dt,dt,ef,ff,gf,res)
 call require(maxval(abs(res))==0.d0,'constant state residual')
 z_residual_scale=1.d0
 ef=0; ff=0; gf=0
 ef(:,2,1,1)=1; ff(:,1,2,1)=2; gf(:,1,1,2)=3
 call calc_R(4,4,nz,1,1,1,dt,dt,dt,ef,ff,gf,res)
 call require(maxval(abs(res/dt-6.d0))<1.d-12,'uniform grid limit')
 print *, 'PASS: stretched mesh, profile, inlet, sponge, extrapolation, metrics and conservation'

contains
 subroutine require(ok,label)
  logical,intent(in)::ok
  character(*),intent(in)::label
  if (.not.ok) then
   print *, 'FAIL: ',label
   stop 1
  endif
 end subroutine
end program
''')
    subprocess.run(['gfortran', '-cpp', '-DOUXSBLI_STRETCHED_Z', '-O0', '-fcheck=all', '-fwrapv', '-fno-range-check', '-ffree-line-length-none',
                    'stubs.f90', 'mod_globals.f90', 'set.f90', 'calc_steps.F90', 'check.f90', '-o', 'check'], cwd=work, check=True)
    subprocess.run([str(work / 'check')], cwd=work, check=True)

    # Also compile the shared routine without the case-specific feature flag.
    subprocess.run(['gfortran', '-cpp', '-ffree-line-length-none', '-c',
                    'calc_steps.F90', '-o', 'legacy_steps.o'], cwd=work, check=True)
