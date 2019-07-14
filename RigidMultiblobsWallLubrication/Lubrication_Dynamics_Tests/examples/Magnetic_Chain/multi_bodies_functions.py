'''
In this module the user can define functions that modified the
code multi_blobs.py. For example, functions to define the
blobs-blobs interactions, the forces and torques on the rigid
bodies or the slip on the blobs.
'''
import numpy as np
import sys
import imp
import os.path
from functools import partial

import utils
from quaternion_integrator.quaternion import Quaternion

# Try to import the forces boost implementation
try:
  import forces_ext
except ImportError:
  pass
# If pycuda is installed import forces_pycuda
try: 
  imp.find_module('pycuda')
  found_pycuda = True
except ImportError:
  found_pycuda = False
if found_pycuda:
  try:
    import pycuda.autoinit
    autoinit_pycuda = True
  except:
    autoinit_pycuda = False
  if autoinit_pycuda:
    import forces_pycuda  

# Override forces_pycuda with user defined functions.
# If forces_pycuda_user_defined does not exists nothing happens.
if found_pycuda:
  forces_pycuda_user_defined = False
  if os.path.isfile('forces_pycuda_user_defined.py'):
    forces_pycuda_user_defined = True
  if forces_pycuda_user_defined:
    del sys.modules['forces_pycuda']
    sys.modules['forces_pycuda'] = __import__('forces_pycuda_user_defined')
    import forces_pycuda
    

def project_to_periodic_image(r, L):
  '''
  Project a vector r to the minimal image representation
  centered around (0,0,0) and of size L=(Lx, Ly, Lz). If 
  any dimension of L is equal or smaller than zero the 
  box is assumed to be infinite in that direction.
  '''
  if L is not None:
    for i in range(3):
      if(L[i] > 0):
        r[i] = r[i] - int(r[i] / L[i] + 0.5 * (int(r[i]>0) - int(r[i]<0))) * L[i]
  return r

def default_zero_r_vectors(r_vectors, *args, **kwargs):
  return np.zeros((r_vectors.size / 3, 3))


def default_zero_blobs(body, *args, **kwargs):
  ''' 
  Return a zero array of shape (body.Nblobs, 3)
  '''
  return np.zeros((body.Nblobs, 3))


def default_zero_bodies(bodies, *args, **kwargs):
  ''' 
  Return a zero array of shape (2*len(bodies), 3)
  '''
  return np.zeros((2*len(bodies), 3))
  

def set_slip_by_ID(body, slip):
  '''
  This function assign a slip function to each body.
  If the body has an associated slip file the function
  "active_body_slip" is assigned (see function below).
  Otherwise the slip is set to zero.

  This function can be override to assign other slip
  functions based on the body ID, (ID of a structure
  is the name of the clones file (without .clones)).
  See the example in
  "examples/pair_active_rods/".
  '''
  if slip is not None:
    active_body_slip_partial = partial(active_body_slip, slip = slip)
    body.function_slip = active_body_slip_partial
  else:
    body.function_slip = default_zero_blobs
  return


def active_body_slip(body, slip):
  '''
  This function set the slip read from the *.slip file to the
  blobs. The slip on the file is given in the body reference 
  configuration (quaternion = (1,0,0,0)) therefore this
  function rotates the slip to the current body orientation.
  
  This function can be used, for example, to model active rods
  that propel along their axes. 
  '''
  # Get rotation matrix
  rotation_matrix = body.orientation.rotation_matrix()

  # Rotate  slip on each blob
  slip_rotated = np.empty((body.Nblobs, 3))
  for i in range(body.Nblobs):
    slip_rotated[i] = np.dot(rotation_matrix, slip[i])
  return slip_rotated


def bodies_external_force_torque(bodies, r_vectors, *args, **kwargs):
  '''
  This function returns the external force-torques acting on the bodies.
  It returns an array with shape (2*len(bodies), 3)
  
  In this is example we just set it to zero.
  '''
  
  ######################
  
  # Change this to add bending and springs and things
  
  ######################
  
  ext_FT = np.zeros((len(bodies), 3))
  blob_radius = kwargs.get('blob_radius')
  eps = kwargs.get('spring_k')
  rest_l = kwargs.get('spring_l')
  k_bend = kwargs.get('k_bend')
  
  r_vecs = [b.location for b in bodies]
  num_parts = len(r_vecs)
  
  # springs
  for k, b in enumerate(bodies):
    if k+1 < num_parts:
      # center spring
      r = r_vecs[k+1] - r_vecs[k]
      r_norm = np.linalg.norm(r)
      force_torque = eps * (r_norm - rest_l) * (r / r_norm) 
      ext_FT[k,:] += force_torque
      ext_FT[k+1,:] -= force_torque
      # twist springs
      #top_k = bodies[k].location + np.dot(bodies[k].orientation.rotation_matrix(),np.array([0.0, 0.0, blob_radius]))
      #top_k1 = bodies[k+1].location + np.dot(bodies[k+1].orientation.rotation_matrix(),np.array([0.0, 0.0, blob_radius]))
      #bot_k = bodies[k].location + np.dot(bodies[k].orientation.rotation_matrix(),np.array([0.0, 0.0, -blob_radius]))
      #bot_k1 = bodies[k+1].location + np.dot(bodies[k+1].orientation.rotation_matrix(),np.array([0.0, 0.0, -blob_radius]))
      
      #top_r = top_k1 - top_k
      #bot_r = bot_k1 - bot_k
      
      #top_r_norm = np.linalg.norm(top_r)
      #force_torque = 0.1*eps * (top_r_norm - rest_l) * (top_r / top_r_norm) 
      #ext_FT[k,:] += force_torque
      #ext_FT[k+1,:] -= force_torque
      
      #bot_r_norm = np.linalg.norm(bot_r)
      #force_torque = 0.1*eps * (bot_r_norm - rest_l) * (bot_r / bot_r_norm) 
      #ext_FT[k,:] += force_torque
      #ext_FT[k+1,:] -= force_torque
      
      
  #bending
  #for k in range(num_parts):
    #if k == 0:
      #Fb = -1.0*k_bend*(-2*r_vecs[1] + r_vecs[0] + r_vecs[2])
    #elif k == 1:
      #Fb = -1.0*k_bend*(-2*r_vecs[0] + 5*r_vecs[1] - 4*r_vecs[2] + r_vecs[3])
    #elif k == (num_parts-1):
      #Fb = -1.0*k_bend*(-2*r_vecs[num_parts-2] + r_vecs[num_parts-1] + r_vecs[num_parts-3])
    #elif k == (num_parts-2):
      #Fb = -1.0*k_bend*(-2*r_vecs[num_parts-1] + 5*r_vecs[num_parts-2] - 4*r_vecs[num_parts-3] + r_vecs[num_parts-4])
    #else:
      #Fb = -1.0*k_bend*(r_vecs[k-2] - 4*r_vecs[k-1] + 6*r_vecs[k] - 4*r_vecs[k+1] + r_vecs[k+2])
    #ext_FT[k,:] += Fb
      
      
  ## bending
  ends = 1.0
  for k in range(num_parts):
    
    k_bend_fn = k_bend*(1.0 - 0.6*(1.0*k/(1.0*num_parts)))
    
    if k == 0:
      vb_0 = (r_vecs[1] - r_vecs[0])
      b_0 = np.linalg.norm(vb_0)
      u_0 = vb_0/b_0 
      
      vb_1 = (r_vecs[2] - r_vecs[1])
      b_1 = np.linalg.norm(vb_1)
      u_1 = vb_1/b_1 
      
      Fb = ends*-1.0*k_bend_fn*(u_1-u_0)/b_0
    elif k == 1:
      vb_0 = (r_vecs[1] - r_vecs[0])
      b_0 = np.linalg.norm(vb_0)
      u_0 = vb_0/b_0 
      
      vb_1 = (r_vecs[2] - r_vecs[1])
      b_1 = np.linalg.norm(vb_1)
      u_1 = vb_1/b_1 
      
      vb_2 = (r_vecs[3] - r_vecs[2])
      b_2 = np.linalg.norm(vb_2)
      u_2 = vb_2/b_2 
      
      Fb = ends*k_bend_fn*( (u_1-u_0)*(1.0/b_0 + 1.0/b_1) - (u_2 - u_1)/b_1 )
    elif k == (num_parts-1):
      vb_N1 = (r_vecs[num_parts-1] - r_vecs[num_parts-2])
      b_N1 = np.linalg.norm(vb_N1)
      u_N1 = vb_N1/b_N1 
      
      vb_N2 = (r_vecs[num_parts-2] - r_vecs[num_parts-3])
      b_N2 = np.linalg.norm(vb_N2)
      u_N2 = vb_N2/b_N2 
      
      Fb = ends*-1.0*k_bend_fn*(u_N1-u_N2)/b_N1
    elif k == (num_parts-2):
      vb_N1 = (r_vecs[num_parts-1] - r_vecs[num_parts-2])
      b_N1 = np.linalg.norm(vb_N1)
      u_N1 = vb_N1/b_N1 
      
      vb_N2 = (r_vecs[num_parts-2] - r_vecs[num_parts-3])
      b_N2 = np.linalg.norm(vb_N2)
      u_N2 = vb_N2/b_N2 
      
      vb_N3 = (r_vecs[num_parts-3] - r_vecs[num_parts-4])
      b_N3 = np.linalg.norm(vb_N3)
      u_N3 = vb_N3/b_N3 
      
      Fb = ends*k_bend_fn*( (u_N1-u_N2)*(1.0/b_N2 + 1.0/b_N1) - (u_N2 - u_N3)/b_N2 )
    else:
      v_k = r_vecs[k+1] - r_vecs[k]
      v_kp1 = r_vecs[k+2] - r_vecs[k+1]
      v_k1 = r_vecs[k] - r_vecs[k-1]
      v_k2 = r_vecs[k-1] - r_vecs[k-2]

      b_k = np.linalg.norm(v_k)
      b_kp1 = np.linalg.norm(v_kp1)
      b_k1 = np.linalg.norm(v_k1)
      b_k2 = np.linalg.norm(v_k2)
      
      u_k = v_k/b_k
      u_kp1 = v_kp1/b_kp1
      u_k1 = v_k1/b_k1
      u_k2 = v_k2/b_k2
      
      
      Fb = k_bend_fn*( (u_k-u_k1)*(1.0/b_k1 + 1.0/b_k) - (u_k1 - u_k2)/b_k1 - (u_kp1 - u_k)/b_k )
    ext_FT[k,:] += Fb
      
      
      
  ## bending
  #for k in range(num_parts):
    #if k == 0:
      #r_ss_1 = (-2*r_vecs[1] + r_vecs[0] + r_vecs[2])/(rest_l**2)
      #Fb = 0*(2.0*k_bend/(rest_l**2))*r_ss_1
    #elif k == 1:
      #r_ss_1 = (-2*r_vecs[1] + r_vecs[0] + r_vecs[2])/(rest_l**2)
      #r_ss_2 = (-2*r_vecs[2] + r_vecs[1] + r_vecs[3])/(rest_l**2)
      #Fb = 0*(k_bend/(rest_l**2))*(-2.0*r_ss_1 + r_ss_2)
    #elif k == (num_parts-1):
      #r_ss_1 = (-2*r_vecs[num_parts-2] + r_vecs[num_parts-1] + r_vecs[num_parts-3])/(rest_l**2)
      #Fb = 0*(2.0*k_bend/(rest_l**2))*r_ss_1
    #elif k == (num_parts-2):
      #r_ss_1 = (-2*r_vecs[num_parts-2] + r_vecs[num_parts-1] + r_vecs[num_parts-3])/(rest_l**2)
      #r_ss_2 = (-2*r_vecs[num_parts-3] + r_vecs[num_parts-2] + r_vecs[num_parts-4])/(rest_l**2)
      #Fb = 0*(k_bend/(rest_l**2))*(-2.0*r_ss_1 + r_ss_2)
    #else:
      #r_ss_L = (-2*r_vecs[k-1] + r_vecs[k-2] + r_vecs[k])/(rest_l**2)
      #r_ss_C = (-2*r_vecs[k] + r_vecs[k-1] + r_vecs[k+1])/(rest_l**2)
      #r_ss_R = (-2*r_vecs[k+1] + r_vecs[k] + r_vecs[k+2])/(rest_l**2)
      #Fb = (k_bend/(rest_l**2))*(-2.0*r_ss_C + r_ss_L + r_ss_R)
    #ext_FT[k,:] += Fb
  
  return ext_FT
  

def blob_external_force(r_vectors, *args, **kwargs):
  '''
  This function compute the external force acting on a
  single blob. It returns an array with shape (3).
  
  In this example we add gravity and a repulsion with the wall;
  the interaction with the wall is derived from the potential

  U(z) = U0 + U0 * (a-z)/b   if z<a
  U(z) = U0 * exp(-(z-a)/b)  iz z>=a

  with 
  e = repulsion_strength_wall
  a = blob_radius
  h = distance to the wall
  b = debye_length_wall
  '''
  f = np.zeros(3)

  # Get parameters from arguments
  blob_mass = kwargs.get('blob_mass')
  blob_radius = kwargs.get('blob_radius')
  g = kwargs.get('g')
  repulsion_strength_wall = kwargs.get('repulsion_strength_wall') 
  debye_length_wall = kwargs.get('debye_length_wall')
  
  repulsion_strength = kwargs.get('repulsion_strength') 
  debye_length = kwargs.get('debye_length')
  
  debye_delta = kwargs.get('debye_delta')
  # Add gravity
  f += -g * blob_mass * np.array([0., 0., 1.0])
  
  # Add wall interaction
  h = r_vectors[2]
  if h > (blob_radius-debye_delta):
    f[2] += (repulsion_strength_wall / debye_length_wall) * np.exp(-(h-(blob_radius-debye_delta))/debye_length_wall)
  else:
    f[2] += (repulsion_strength_wall / debye_length_wall)
    
  # Add wall interaction
  if h > (blob_radius):
    f[2] += (repulsion_strength / debye_length) * np.exp(-(h-blob_radius)/debye_length)
  else:
    f[2] += (repulsion_strength / debye_length)
      
  return f


def calc_one_blob_forces(r_vectors, *args, **kwargs):
  '''
  Compute one-blob forces. It returns an array with shape (Nblobs, 3).
  '''
  Nblobs = r_vectors.size / 3
  force_blobs = np.zeros((Nblobs, 3))
  r_vectors = np.reshape(r_vectors, (Nblobs, 3))
  
  # Loop over blobs
  for blob in range(Nblobs):
    force_blobs[blob] += blob_external_force(r_vectors[blob], *args, **kwargs)   

  return force_blobs


def set_blob_blob_forces(implementation):
  '''
  Set the function to compute the blob-blob forces
  to the right function.
  The implementation in pycuda is much faster than the
  one in C++, which is much faster than the one python; 
  To use the pycuda implementation is necessary to have 
  installed pycuda and a GPU with CUDA capabilities. To
  use the C++ implementation the user has to compile 
  the file blob_blob_forces_ext.cc.   
  '''
  if implementation == 'None':
    return default_zero_r_vectors
  elif implementation == 'python':
    return calc_blob_blob_forces_python
  elif implementation == 'C++':
    return calc_blob_blob_forces_boost 
  elif implementation == 'pycuda':
    return forces_pycuda.calc_blob_blob_forces_pycuda


def blob_blob_force(r, *args, **kwargs):
  '''
  This function compute the force between two blobs
  with vector between blob centers r.

  In this example the force is derived from the potential
  
  U(r) = U0 + U0 * (2*a-r)/b   if z<2*a
  U(r) = U0 * exp(-(r-2*a)/b)  iz z>=2*a
  
  with
  eps = potential strength
  r_norm = distance between blobs
  b = Debye length
  a = blob_radius
  '''
  # Get parameters from arguments
  L = kwargs.get('periodic_length')
  a = kwargs.get('blob_radius')
  time_s = kwargs.get('time_s')
  
  repulsion_strength_wall = kwargs.get('repulsion_strength_wall') 
  debye_length_wall = kwargs.get('debye_length_wall')
  
  repulsion_strength = kwargs.get('repulsion_strength') 
  debye_length = kwargs.get('debye_length')
  
  debye_delta = kwargs.get('debye_delta')

  # Compute force
  project_to_periodic_image(r, L)
  r_norm = np.linalg.norm(r)
  r_hat = r/r_norm
  
  ################
  #phi = np.arctan2(r_hat[1],r_hat[2])
  #theta = np.arccos(r_hat[0])
  
  #theta_hat = np.array([-1.0*np.sin(theta), np.cos(theta)*np.sin(phi), np.cos(theta)*np.cos(phi)])
  
  ################
  ##C = kwargs.get('mag_force')
  ##force_torque = C*((2*a/r_norm)**4)*( (3*(np.cos(theta))**2 - 1.0)*r_hat + (np.sin(2*theta))*theta_hat )
  ################
  
  ################
  #C = kwargs.get('mag_force')
  #Bratio = 2.94/1.43 # = Bxy/Bz
  #beta = 0.25
  #chi = beta*(a/r_norm)**3
  
  #### angles ###
  #phi_z = np.arctan2(1.5*chi*np.sin(2.0*theta), (1 + 0.5*chi* (3.0*np.cos(2.0*theta) - 1.0) ) )
  #phi_xy = np.arctan2(1.5*chi*np.sin(2.0*theta), (1 - 0.5*chi* (3.0*np.cos(2.0*theta) + 1.0) ) )
  ###############
  #Fz_C = -C*((a/r_norm)**4)
  #Fz_chi = ( (1 - chi + (5.0/2.0)*chi**2 + 3*chi*(1 - 0.5*chi)*np.cos(2*theta)) / ( ((1 + chi)**2) * ((1 - 2*chi)**2) ) )
  #Fz = Fz_C * Fz_chi * ( ((3*np.cos(2 * (theta - phi_z)) + 1)/2.0) * r_hat + (np.sin(2*(theta-phi_z))) * theta_hat )
	
  
  
  
  #Fxy_C = 0.5*(Bratio**2)*C*((a/r_norm)**4)
  #Fxy_chi = ( (1 - chi + (5.0/2.0)*chi**2 - 3*chi*(1 - 0.5*chi)*np.cos(2*theta)) / ( ((1 + chi)**2) * ((1 - 2*chi)**2) ) )
  #Fxy = Fxy_C * Fxy_chi * ( ((1 - 3*np.cos(2 * (theta - phi_xy)) + 1)/2.0) * r_hat - (np.sin(2*(theta-phi_xy))) * theta_hat ) - Fxy_C * (1.0/(1 + chi)**2) * r_hat
  
  #force_torque = Fz + Fxy
  ###############
  
  #################
  C = kwargs.get('mag_force')
  Bratio = 2.94/1.43 # = Bxy/Bz
  #Mom = np.array([1.43, 2.94*np.sin(2*np.pi*50*time_s), 2.94*np.cos(2*np.pi*50*time_s)])
  #Mom = np.array([1.3, 2.2*np.sin(2*np.pi*50*time_s), 2.2*np.cos(2*np.pi*50*time_s)])
  Mom = np.array([0.0, 2.2*np.sin(2*np.pi*50*time_s), 2.2*np.cos(2*np.pi*50*time_s)])
  m0 = np.linalg.norm(Mom)
  m = Mom/m0
  
  force_torque = 0.0*r_hat
  if r_norm < 10*a:
    force_torque += C * ((2.0*a/r_norm)**4) * ( 2.0 * (np.dot(m,r_hat)) * m - (5.0 * ((np.dot(m,r_hat))**2) - 1.0) * r_hat )
  #################
  
  
  #print "no interparticle force"
  offset = 2.0*a-2.0*debye_delta
  if r_norm > (offset):
    force_torque += -((repulsion_strength_wall / debye_length_wall) * np.exp(-(r_norm-(offset)) / debye_length_wall) / np.maximum(r_norm, np.finfo(float).eps)) * r 
  else:
    force_torque += -((repulsion_strength_wall / debye_length_wall) / np.maximum(r_norm, np.finfo(float).eps)) * r;
    
  offset = 2.0*a
  if r_norm > (offset):
    force_torque += -((repulsion_strength / debye_length) * np.exp(-(r_norm-(offset)) / debye_length) / np.maximum(r_norm, np.finfo(float).eps)) * r 
  else:
    force_torque += -((repulsion_strength / debye_length) / np.maximum(r_norm, np.finfo(float).eps)) * r;
  
  
  return force_torque
  

def calc_blob_blob_forces_python(r_vectors, *args, **kwargs):
  '''
  This function computes the blob-blob forces and returns
  an array with shape (Nblobs, 3).
  '''
  Nblobs = r_vectors.size / 3
  force_blobs = np.zeros((Nblobs, 3))

  # Double loop over blobs to compute forces
  for i in range(Nblobs-1):
    for j in range(i+1, Nblobs):
      # Compute vector from j to u
      r = r_vectors[j] - r_vectors[i]
      force = blob_blob_force(r, *args, **kwargs)
      force_blobs[i] += force
      force_blobs[j] -= force

  return force_blobs


def calc_blob_blob_forces_boost(r_vectors, *args, **kwargs):
  '''
  Call a boost function to compute the blob-blob forces.
  '''
  # Get parameters from arguments
  L = kwargs.get('periodic_length')
  eps = kwargs.get('repulsion_strength')
  b = kwargs.get('debye_length')  
  blob_radius = kwargs.get('blob_radius')  

  number_of_blobs = r_vectors.size / 3
  r_vectors = np.reshape(r_vectors, (number_of_blobs, 3))
  forces = np.empty(r_vectors.size)
  if L is None:
    L = -1.0*np.ones(3)

  forces_ext.calc_blob_blob_forces(r_vectors, forces, eps, b, blob_radius, number_of_blobs, L)
  return np.reshape(forces, (number_of_blobs, 3))


def set_body_body_forces_torques(implementation):
  '''
  Set the function to compute the body-body forces
  to the right function. 
  '''
  if implementation == 'None':
    return default_zero_bodies
  elif implementation == 'python':
    return calc_body_body_forces_torques_python


def body_body_force_torque(r, quaternion_i, quaternion_j, *args, **kwargs):
  '''
  This function compute the force between two bodies
  with vector between locations r.
  In this example the torque is zero and the force 
  is derived from a Yukawa potential
  
  U = eps * exp(-r_norm / b) / r_norm
  
  with
  eps = potential strength
  r_norm = distance between bodies' location
  b = Debye length
  '''
  force_torque = np.zeros((2, 3))

  # Get parameters from arguments
  L = kwargs.get('periodic_length')
  eps = kwargs.get('repulsion_strength')
  b = kwargs.get('debye_length')
  
  return force_torque


def calc_body_body_forces_torques_python(bodies, r_vectors, *args, **kwargs):
  '''
  This function computes the body-body forces and torques and returns
  an array with shape (2*Nblobs, 3).
  '''
  Nbodies = len(bodies)
  force_torque_bodies = np.zeros((2*len(bodies), 3))
  
  # Double loop over bodies to compute forces
  for i in range(Nbodies-1):
    for j in range(i+1, Nbodies):
      # Compute vector from j to u
      r = bodies[j].location - bodies[i].location
      force_torque = body_body_force_torque(r, bodies[i].orientation, bodies[j].orientation, *args, **kwargs)
      # Add forces
      force_torque_bodies[2*i] += force_torque[0]
      force_torque_bodies[2*j] -= force_torque[0]
      # Add torques
      force_torque_bodies[2*i+1] += force_torque[1]
      force_torque_bodies[2*j+1] -= force_torque[1]

  return force_torque_bodies


def force_torque_calculator_sort_by_bodies(bodies, r_vectors, *args, **kwargs):
  '''
  Return the forces and torque in each body with
  format [f_1, t_1, f_2, t_2, ...] and shape (2*Nbodies, 3),
  where f_i and t_i are the force and torque on the body i.
  '''
  # Create auxiliar variables
  Nblobs = r_vectors.size / 3
  force_torque_bodies = np.zeros((2*len(bodies), 3))
  force_blobs = np.zeros((Nblobs, 3))
  blob_mass = 1.0
  blob_radius = bodies[0].blob_radius

  # Compute one-blob forces (same function for all blobs)
  force_blobs += calc_one_blob_forces(r_vectors, blob_radius = blob_radius, blob_mass = blob_mass, *args, **kwargs)

  # Compute blob-blob forces (same function for all pair of blobs)
  force_blobs += calc_blob_blob_forces(r_vectors, blob_radius = blob_radius, *args, **kwargs)  

  force_blobs += bodies_external_force_torque(bodies, r_vectors, blob_radius = blob_radius, *args, **kwargs)

  # Compute body force-torque forces from blob forces
  offset = 0
  for k, b in enumerate(bodies):
    # Add force to the body
    force_torque_bodies[2*k:(2*k+1)] += sum(force_blobs[offset:(offset+b.Nblobs)])
    # Add torque to the body
    R = b.calc_rot_matrix()  
    force_torque_bodies[2*k+1:2*k+2] += np.dot(R.T, np.reshape(force_blobs[offset:(offset+b.Nblobs)], 3*b.Nblobs))
    offset += b.Nblobs

  # Add one-body external force-torque
  #force_torque_bodies += bodies_external_force_torque(bodies, r_vectors, blob_radius = blob_radius, *args, **kwargs)

  # Add body-body forces (same for all pair of bodies)
  force_torque_bodies += calc_body_body_forces_torques(bodies, r_vectors, *args, **kwargs)
  return force_torque_bodies


def preprocess(bodies, *args, **kwargs):
  '''
  This function is call at the start of the schemes.
  The default version do nothing, it should be modify by
  the user if he wants to change the schemes.
  '''
  return

def postprocess(bodies, *args, **kwargs):
  '''
  This function is call at the end of the schemes but
  before checking if the postions are a valid configuration.
  The default version do nothing, it should be modify by
  the user if he wants to change the schemes.
  '''
  return


# Override force interactions by user defined functions.
# This only override the functions implemented in python.
# If user_defined_functions is empty or does not exists
# this import does nothing.
user_defined_functions_found = False
if os.path.isfile('user_defined_functions.py'):
  user_defined_functions_found = True
if user_defined_functions_found:
  import user_defined_functions

