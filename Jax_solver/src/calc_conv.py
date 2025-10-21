import jax
import jax.numpy as jnp


@jax.jit
def calc_E(gamma:float, rho:jnp.float64, u:jnp.float64, v:jnp.float64, w:jnp.float64, p:jnp.float64):
  # rho, u, v, w, p [nz-2,ny-2,nx]
  E1 = 0.25e0 * (rho[:,:,:-1] + rho[:,:,1:]) * (u[:,:,:-1] + u[:,:,1:])
  E2 = E1 * 0.5e0 * (u[:,:,:-1] + u[:,:,1:]) + 0.5e0 * (p[:,:,:-1] + p[:,:,1:])
  E3 = E1 * 0.5e0 * (v[:,:,:-1] + v[:,:,1:])
  E4 = E1 * 0.5e0 * (w[:,:,:-1] + w[:,:,1:])
  E5 = E1 * 0.5e0 * (p[:,:,:-1] / rho[:,:,:-1] + p[:,:,1:] / rho[:,:,1:]) / (gamma - 1.e0) \
          + 0.5e0 * (u[:,:,:-1] * p[:,:,1:] + u[:,:,1:] * p[:,:,:-1]) \
          + 0.5e0 * E1 * (u[:,:,:-1] * u[:,:,1:] + v[:,:,:-1] * v[:,:,1:] + w[:,:,:-1] * w[:,:,1:])
  return jnp.stack([E1, E2, E3, E4, E5], axis=-1)


@jax.jit
def calc_F(gamma:float, rho:jnp.float64, u:jnp.float64, v:jnp.float64, w:jnp.float64, p:jnp.float64):
  # rho, u, v, w, p [nz-2,ny,nx-2]
  F1 = 0.25e0 * (rho[:,:-1,:] + rho[:,1:,:]) * (v[:,:-1,:] + v[:,1:,:])
  F2 = F1 * 0.5e0 * (u[:,:-1,:] + u[:,1:,:])
  F3 = F1 * 0.5e0 * (v[:,:-1,:] + v[:,1:,:]) + 0.5e0 * (p[:,:-1,:] + p[:,1:,:])
  F4 = F1 * 0.5e0 * (w[:,:-1,:] + w[:,1:,:])
  F5 = F1 * 0.5e0 * (p[:,:-1,:] / rho[:,:-1,:] + p[:,1:,:] / rho[:,1:,:]) / (gamma - 1.e0) \
          + 0.5e0 * (v[:,:-1,:] * p[:,1:,:] + v[:,1:,:] * p[:,:-1,:]) \
          + 0.5e0 * F1 * (u[:,:-1,:] * u[:,1:,:] + v[:,:-1,:] * v[:,1:,:] + w[:,:-1,:] * w[:,1:,:])
  return jnp.stack([F1, F2, F3, F4, F5], axis=-1)


@jax.jit
def calc_G(gamma:float, rho:jnp.float64, u:jnp.float64, v:jnp.float64, w:jnp.float64, p:jnp.float64):
  # rho, u, v, w, p [nz,ny-2,nx-2]
  G1 = 0.25e0 * (rho[:-1,:,:] + rho[1:,:,:]) * (w[:-1,:,:] + w[1:,:,:])
  G2 = G1 * 0.5e0 * (u[:-1,:,:] + u[1:,:,:])
  G3 = G1 * 0.5e0 * (v[:-1,:,:] + v[1:,:,:])
  G4 = G1 * 0.5e0 * (w[:-1,:,:] + w[1:,:,:]) + 0.5e0 * (p[:-1,:,:] + p[1:,:,:])
  G5 = G1 * 0.5e0 * (p[:-1,:,:] / rho[:-1,:,:] + p[1:,:,:] / rho[1:,:,:]) / (gamma - 1.e0) \
          + 0.5e0 * (w[:-1,:,:] * p[1:,:,:] + w[1:,:,:] * p[:-1,:,:]) \
          + 0.5e0 * G1 * (u[:-1,:,:] * u[1:,:,:] + v[:-1,:,:] * v[1:,:,:] + w[:-1,:,:] * w[1:,:,:])
  return jnp.stack([G1, G2, G3, G4, G5], axis=-1)

