import argparse
import numpy as np
import scipy.linalg
import subprocess
import cPickle
from functools import partial
import sys
import time
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
  
  #h=1.260716491714539
  #L_wall = Lub.Lub_Resist_Wall(h, eta, a)
  #Res_Wall_Coeffs = Lub.Test_Single_Blob_Resistance(eta, a, 500, 0.0009998)
  
  #Res_Pair_Coeffs = Lub.Pair_Blob_Resistance_Free_Space(eta, a, 3001, 0.0009998)
  Res_Pair_Coeffs = Lub.Pair_Blob_Resistance_Free_Space_Test(eta, a, 3001, 0.0009998)
  
  #LubNorms = []
  #print len(bodies)
  #for k in range(0,len(bodies),2):
    #norm = np.linalg.norm(Lub.Resist(bodies[k].location, bodies[k+1].location, eta, a))
    #LubNorms.append(norm)
   
  #print len(LubNorms)
  #for i, ID in enumerate(structures_ID):
    #name = 'Lub_Sup_test' + '.' + ID + '.errors'
    #with open(name, 'w') as f_ID:
      #f_ID.write('%s \n' % (LubNorms[0]))
    #status = 'a'
    #with open(name, status) as f_ID:
      #for j in range(1,len(LubNorms)):
        #f_ID.write('%s \n' % (LubNorms[j]))

  ##start = time.time()
  ##pair_lub = Lub.pair_resistance_sup_tree(eta, a, 2.5)
  ##end = time.time()
  ##print 'time for pair resistance cKD: '+ str((end - start))
  ##start = time.time()
  ##pair_lub_llc = Lub.pair_resistance_sup_llc(eta, a, 2.5)
  ##end = time.time()
  ##print 'time for pair resistance llc: '+ str((end - start))
  
  ##F = np.random.rand(6*num_bodies)
  
  
  ##start = time.time()
  ##pair_lub_m = Lub.pair_resistance_mult_tree(F, eta, a, 2.5)
  ##end = time.time()
  ##print 'time for cKD mult: ' + str((end - start)) + '. Result:  ' + str(np.linalg.norm(pair_lub_m))
  ##start = time.time()
  ##pair_lub = Lub.pair_resistance_sup_tree(eta, a, 2.5)
  ##pair_lub_mm = np.dot(pair_lub,F)
  ##end = time.time()
  ##print 'time for mat mult: '+ str((end - start)) + '. Result:  ' + str(np.linalg.norm(pair_lub_mm))
  
  
  ##start = time.time()
  ##F_det = 1.0*np.arange(6*num_bodies)
  ##pair_lub_mb = Lub.pair_resistance_mult_MB_tree(F_det, eta, a, 2.5)
  ##end = time.time()
  ##print 'time for MB mult: '+ str((end - start)) + '. Result:  ' + str(np.linalg.norm(pair_lub_mb))
  
  
  ##r_locs = np.array([b.location for b in bodies])
  ##print 'min body height = ' + str(np.min(r_locs[:,2]))
  
  ##print np.linalg.norm(pair_lub)
  ##print np.linalg.norm(pair_lub_llc)
  ##name = 'Lub_Sup_MATLAB_Compare.dat'
  ##np.savetxt(name, pair_lub, delimiter='  ')
  
  #new_bodies = []
  #new_bodies.append(bodies[0])
  #new_bodies.append(bodies[1])
  #new_bodies[0].location = np.array([0.0, 0.0, 5.0])
  #new_bodies[1].location = np.array([2.02, 0.0, 5.0])

  #print(new_bodies[0].location)
  #print(new_bodies[1].location)
  
  #Lub.bodies = new_bodies
  #num_bodies = 2
  
  
  #########F_1 = np.random.rand(6*num_bodies)
  #########F_2 = np.random.rand(6*num_bodies)
  #########F = np.concatenate((0.0*F_1,F_2))
  #########start = time.time()
  #########UW = Lub.Lubrucation_operator_mult(F, eta, a, 2.5)
  #########end = time.time()
  #########print 'time for operator apply: '+ str((end - start)) + '. Result:  ' + str(np.linalg.norm(UW))
  
  #########(U_alt, alt_res) = Lub.Lubrucation_alt_solve(-1.0*F_2, eta, a, 2.5)
  #########(U_gmres, saddle_res) = Lub.Lubrucation_solve(F, eta, a, 2.5)
  
  #########plt.hold(True)
  #########plt.semilogy(alt_res,'-ro')
  #########plt.semilogy(saddle_res,'-bo')
  #########plt.show()
  
  #########print(np.linalg.norm(U_gmres[6*num_bodies::]-U_alt)/np.linalg.norm(U_alt))
  
  #np.linalg.norm(Lub.Resist(bodies[0].location, bodies[1].location, eta, a,inv=True))
  #print np.linalg.norm(Lub.Resist(bodies[0].location, bodies[1].location, eta, a))
  
  #s1 = np.array([0.0, 0.0, 0.0])
  #s2 = np.random.rand(3)
  #s2 = 2.101*(s2/np.linalg.norm(s2))
  
  #start = time.time()
  #print np.linalg.norm(Lub.Resist_Method(s1, s2, eta, a,'JO'))
  #end = time.time()
  #print(end - start)
  #start = time.time()
  #print np.linalg.norm(Lub.Resist_Method(s1, s2, eta, a,'JOR'))
  #end = time.time()
  #print(end - start)
  
  