'''
Use this module to override forces interactions defined in
multi_body_functions.py. See an example in the file
RigidMultiblobsWall/multi_bodies/examples/user_defined_functions.py



In this module we override the default blob-blob, blob-wall and
body-body interactions used by the code. To use this implementation
copy this file to
RigidMultiblobsWall/multi_bodies/user_defined_functions.py


This module defines (and override) the following interactions:

1. blob-wall forces, they are derived from the potential:
  U = e * a * exp(-(h-a) / b) / (h - a)
  with
  e = repulsion_strength_wall
  a = blob_radius
  h = distance to the wall
  b = debye_length_wall

2. blob-blob forces, they are derived from the potential:
  U = eps * exp(-r_norm / b) / r_norm
  with
  eps = potential strength
  r_norm = distance between blobs
  b = Debye length

3. body-body forces and torques. The torques are zero and the forces
  are derived from the potential:
  U = 0.5 * eps * (r_norm - 1.0)**2
  with
  eps = potential strength
  r_norm = distance between bodies' location
'''

import multi_bodies_functions
from multi_bodies_functions import *
import math as m



# Override blob_external_force
def blob_external_force_new(r_vectors, *args, **kwargs):
  '''
  This function compute the external force acting on a
  single blob. It returns an array with shape (3).

  In this example we add gravity and a repulsion with the wall;
  the interaction with the wall is derived from a Yukawa-like
  potential
  U = e * a * exp(-(h-a) / b) / (h - a)
  with
  e = repulsion_strength_wall
  a = blob_radius
  h = distance to the wall
  b = debye_length_wall
  '''
  f = np.zeros(3)

  return f


def bodies_external_force_torque_new(bodies, r_vectors, *args, **kwargs):
  '''
  This function returns the external force-torques acting on the bodies.
  It returns an array with shape (2*len(bodies), 3)
  
  In this is example we just set it to zero.
  '''
  part = kwargs.get('part')
  comp = kwargs.get('comp')
  sign = kwargs.get('sign')
  
  
  f = np.zeros((2*len(bodies), 3))

  f[2*part,comp] = sign*1.0

  return f

multi_bodies_functions.blob_external_force = blob_external_force_new
multi_bodies_functions.bodies_external_force_torque = bodies_external_force_torque_new