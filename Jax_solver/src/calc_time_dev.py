import jax
import jax.numpy as jnp
import sys
sys.path.append("../KV")
from set import set_bc
from calc_conv import calc_E, calc_F, calc_G
from calc_visc import calc_Ev, calc_Fv, calc_Gv


@jax.jit
def calcR(dt:float, dx:jnp.float64, dy:jnp.float64, dz:jnp.float64, \
          E:jnp.float64, F:jnp.float64, G:jnp.float64):
  R = dt * (dy[None,:-1,None,None] * dz[:-1,None,None,None] * (-E[:,:,:-1,:] + E[:,:,1:,:]) \
          + dz[:-1,None,None,None] * dx[None,None,:-1,None] * (-F[:,:-1,:,:] + F[:,1:,:,:]) \
          + dx[None,None,:-1,None] * dy[None,:-1,None,None] * (-G[:-1,:,:,:] + G[1:,:,:,:]))
  return R


@jax.jit
def step1(dt:float, dx:jnp.float64, dy:jnp.float64, dz:jnp.float64, \
          E:jnp.float64, F:jnp.float64, G:jnp.float64, Q:jnp.float64):
  R  = calcR(dt, dx, dy, dz, E, F, G)
  Q2 = jnp.zeros_like(Q)
  Q2 = Q2.at[1:-1,1:-1,1:-1,:].set(Q[1:-1,1:-1,1:-1,:] - R)
  return Q2


@jax.jit
def step2(dt:float, dx:jnp.float64, dy:jnp.float64, dz:jnp.float64, \
          E:jnp.float64, F:jnp.float64, G:jnp.float64, Q:jnp.float64, Q2:jnp.float64):
  R  = calcR(dt, dx, dy, dz, E, F, G)
  Q3 = jnp.zeros_like(Q)
  Q3 = Q3.at[1:-1,1:-1,1:-1,:].set((0.25e0 * (3.e0 * Q[1:-1,1:-1,1:-1,:] + Q2[1:-1,1:-1,1:-1,:] - R)))
  return Q3


@jax.jit
def step3(dt:float, dx:jnp.float64, dy:jnp.float64, dz:jnp.float64, \
          E:jnp.float64, F:jnp.float64, G:jnp.float64, Q:jnp.float64, Q3:jnp.float64):
  R = calcR(dt, dx, dy, dz, E, F, G)
  Q = Q.at[1:-1,1:-1,1:-1,:].set((Q[1:-1,1:-1,1:-1,:] + 2.e0 * Q3[1:-1,1:-1,1:-1,:] - 2.e0 * R) / 3.e0)
  return Q


@jax.jit
def calc_EFG(dx:jnp.float64, dy:jnp.float64, dz:jnp.float64, gamma:float, Rgas:float, Cp:float, Pr:float, J:jnp.float64, Q:jnp.float64):
  rho = Q[:,:,:,0] * J
  u   = Q[:,:,:,1] / Q[:,:,:,0]
  v   = Q[:,:,:,2] / Q[:,:,:,0]
  w   = Q[:,:,:,3] / Q[:,:,:,0]
  p   = (gamma - 1.e0) * (Q[:,:,:,4] * J  - 0.5e0 * rho * (u**2 + v**2 + w**2))

  E = calc_E(gamma, rho[1:-1,1:-1,:], u[1:-1,1:-1,:], v[1:-1,1:-1,:], w[1:-1,1:-1,:], p[1:-1,1:-1,:])
  F = calc_F(gamma, rho[1:-1,:,1:-1], u[1:-1,:,1:-1], v[1:-1,:,1:-1], w[1:-1,:,1:-1], p[1:-1,:,1:-1])
  G = calc_G(gamma, rho[:,1:-1,1:-1], u[:,1:-1,1:-1], v[:,1:-1,1:-1], w[:,1:-1,1:-1], p[:,1:-1,1:-1])
  E = calc_Ev(dx, dy, dz, rho, u, v, w, p, Rgas, Cp, Pr, E)
  F = calc_Fv(dx, dy, dz, rho, u, v, w, p, Rgas, Cp, Pr, F)
  G = calc_Gv(dx, dy, dz, rho, u, v, w, p, Rgas, Cp, Pr, G)
  return E, F, G


def Runge_Kutta(itr, x):
  gamma, Rgas, Cp, Pr, M0, rho0, u0, p0, T0, dt, dx, dy, dz, J, Q = x
  
  # 1st step
  E, F, G = calc_EFG(dx, dy, dz, gamma, Rgas, Cp, Pr, J, Q)
  Qs = step1(dt, dx, dy, dz, E, F, G, Q)
  Qs = set_bc(gamma, Rgas, M0, rho0, u0, p0, T0, J, Qs)
  
  # 2nd step
  E, F, G = calc_EFG(dx, dy, dz, gamma, Rgas, Cp, Pr, J, Qs)
  Qs = step2(dt, dx, dy, dz, E, F, G, Q, Qs)
  Qs = set_bc(gamma, Rgas, M0, rho0, u0, p0, T0, J, Qs)
  
  # 3rd step
  E, F, G = calc_EFG(dx, dy, dz, gamma, Rgas, Cp, Pr, J, Qs)
  Q = step3(dt, dx, dy, dz, E, F, G, Q, Qs)
  Q = set_bc(gamma, Rgas, M0, rho0, u0, p0, T0, J, Q)
  return (gamma, Rgas, Cp, Pr, M0, rho0, u0, p0, T0, dt, dx, dy, dz, J, Q)

