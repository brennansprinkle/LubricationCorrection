import argparse
import numpy as np
import scipy.linalg
import subprocess
import cPickle
from functools import partial
import sys
import time
import copy
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

# Find project functions
found_functions = False
path_to_append = ''
sys.path.append('../Lubrication')
while found_functions is False:
  try:
    from Lubrication import Lubrication
    from body import body 
    from read_input import read_input
    from read_input import read_vertex_file
    from read_input import read_clones_file
    from read_input import read_slip_file
    import utils
    found_functions = True
  except ImportError:
    path_to_append += '../'
    print 'searching functions in path ', path_to_append
    sys.path.append(path_to_append)
    if len(path_to_append) > 21:
      print '\nProjected functions not found. Edit path in multi_bodies.py'
      sys.exit()

if __name__ == '__main__':
  # Get command line arguments
  parser = argparse.ArgumentParser(description='Run a multi-body simulation and save trajectory.')
  parser.add_argument('--input-file', dest='input_file', type=str, default='data.main', help='name of the input file')
  parser.add_argument('--print-residual', action='store_true', help='print gmres and lanczos residuals')
  args=parser.parse_args()
  input_file = args.input_file

  # Read input file
  read = read_input.ReadInput(input_file)
   
  # Set some variables for the simulation
  eta = read.eta 
  a = read.blob_radius
  output_name = read.output_name 
  structures = read.structures
  structures_ID = read.structures_ID

  # Copy input file to output
  subprocess.call(["cp", input_file, output_name + '.inputfile'])

  # Set random generator state
  if read.random_state is not None:
    with open(read.random_state, 'rb') as f:
      np.random.set_state(cPickle.load(f))
  elif read.seed is not None:
    np.random.seed(int(read.seed))
  
  # Save random generator state
  with open(output_name + '.random_state', 'wb') as f:
    cPickle.dump(np.random.get_state(), f)

  # Create rigid bodies
  bodies = []
  body_types = []
  body_names = []
  for ID, structure in enumerate(structures):
    print 'Creating structures = ', structure[1]
    # Read vertex and clones files
    struct_ref_config = read_vertex_file.read_vertex_file(structure[0])
    num_bodies_struct, struct_locations, struct_orientations = read_clones_file.read_clones_file(structure[1])
    # Read slip file if it exists
    slip = None
    if(len(structure) > 2):
      slip = read_slip_file.read_slip_file(structure[2])
    body_types.append(num_bodies_struct)
    body_names.append(structures_ID[ID])
    # Create each body of type structure
    for i in range(num_bodies_struct):
      b = body.Body(struct_locations[i], struct_orientations[i], struct_ref_config, a)
      b.ID = structures_ID[ID]
      # Calculate body length for the RFD
      if i == 0:
        b.calc_body_length()
      else:
        b.body_length = bodies[-1].body_length
      # Append bodies to total bodies list
      bodies.append(b)
  bodies = np.array(bodies)

  # Set some more variables
  num_of_body_types = len(body_types)
  num_bodies = bodies.size
  Nblobs = sum([x.Nblobs for x in bodies])
  
  Lub=Lubrication(bodies)
  Lub.load_WS_coefficient_interp_data()
  Lub.set_WS_coefficient_interp_functions() # use kind = 'cubic' for WS to get closer to the linear resistance result
  Lub.load_JO_coefficient_interp_data()
  Lub.set_JO_coefficient_interp_functions() # use kind = 'cubic' for JO to get closer to the series result
  Lub.load_MB_wall_interp_data()
  Lub.set_MB_wall_interp_functions()
  Lub.load_MB_coefficient_interp_data(1)
  Lub.set_MB_coefficient_interp_functions()
  Lub.load_res_MB_wall_interp_data()
  Lub.set_MB_res_wall_interp_functions()
  
  
  for dist in [0.025,0.25,1.0,4.0,6.0]:
    d_count = -1
    for h in range(51):
      height = 1.0 + (h+1.0)*0.01
      d_count += 1
      new_bodies = []
      for k, b in enumerate(bodies):
	n_b = copy.deepcopy(b)
	new_bodies.append(n_b)
	new_bodies[k].location[2] -= 2.0
	factor = (2*a + dist*a)/3.0
	new_bodies[k].location *= factor
	new_bodies[k].location[2] += height
	print(new_bodies[k].location)
      Lub.bodies = new_bodies
      F = np.zeros(4*6)
      F[14] = 1.0;
      UW, res = Lub.Lubrucation_alt_solve(F, eta, a, 4.5)
      if d_count == 0:
	f = open("tetra_particle_3_perp_velocities_Lub_dist_" + str(dist) + ".dat", "w")
	np.savetxt(f, np.array([height]), newline=" ")
	f = open("tetra_particle_3_perp_velocities_Lub_dist_" + str(dist) + ".dat", "a")
	np.savetxt(f, UW, newline=" ")
	f.write('\n')
      else:
	f = open("tetra_particle_3_perp_velocities_Lub_dist_" + str(dist) + ".dat", "a")
	np.savetxt(f, np.array([height]), newline=" ")
	np.savetxt(f, UW, newline=" ")
	f.write('\n')

  
  