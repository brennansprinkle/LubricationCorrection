import argparse
import numpy as np
import scipy.linalg
import scipy.sparse.linalg as spla
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
        print
        'searching functions in path ', path_to_append
        sys.path.append(path_to_append)
        if len(path_to_append) > 21:
            print
            '\nProjected functions not found. Edit path in multi_bodies.py'
            sys.exit()

if __name__ == '__main__':
    # Get command line arguments
    parser = argparse.ArgumentParser(description='Run a multi-body simulation and save trajectory.')
    parser.add_argument('--input-file', dest='input_file', type=str, default='data.main', help='name of the input file')
    parser.add_argument('--print-residual', action='store_true', help='print gmres and lanczos residuals')
    args = parser.parse_args()
    input_file = args.input_file

    # Read input file
    read = read_input.ReadInput(input_file)

    # Set some variables for the simulation
    eta = read.eta
    a = read.blob_radius
    output_name = read.output_name
    structures = read.structures
    print structures
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

    for cuts in [100,200,400,1600,3200]: #range(225,2500,225) range(100,1100,100) range(400,4400,400)
      bfname = 'run_phi_0.4_n_'+str(cuts) #'tri_lattice_rad_squared_'+str(cuts)+'_dist_6'
      # Create rigid bodies
      bodies = []
      body_types = []
      body_names = []
      for ID, structure in enumerate(structures):
	  structure[1] = 'Test_Configs_2/'+bfname+'.clones'
	  print 'Creating structures = ', structure[1]
	  # Read vertex and clones files
	  struct_ref_config = read_vertex_file.read_vertex_file(structure[0])
	  num_bodies_struct, struct_locations, struct_orientations = read_clones_file.read_clones_file(structure[1])
	  # Read slip file if it exists
	  slip = None
	  if (len(structure) > 2):
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

      Lub = Lubrication(bodies)
      Lub.load_WS_coefficient_interp_data()
      Lub.set_WS_coefficient_interp_functions()  # use kind = 'cubic' for WS to get closer to the linear resistance result
      Lub.load_JO_coefficient_interp_data()
      Lub.set_JO_coefficient_interp_functions()  # use kind = 'cubic' for JO to get closer to the series result
      Lub.load_MB_wall_interp_data()
      Lub.set_MB_wall_interp_functions()
      Lub.load_MB_coefficient_interp_data(1)
      Lub.set_MB_coefficient_interp_functions()
      Lub.load_res_MB_wall_interp_data()
      Lub.set_MB_res_wall_interp_functions()


      #linear_operator_partial = partial(Lub.M_DR_Mult, eta=eta, a=a, cutoff=4.5)
      #A = spla.LinearOperator((6*num_bodies, 6*num_bodies), matvec = linear_operator_partial, dtype='float64')
      #print 'start'
      #start = time.time()
      #vals, vecs = spla.eigsh(A, k=20,tol=1e-5,which='LM')
      #end = time.time()
      #print 'time for eigs: '+ str((end - start))
      #print vals
      
      #linear_operator_partial = partial(Lub.M_DR_Mult, eta=eta, a=a, cutoff=4.5)
      #A = spla.LinearOperator((6*num_bodies, 6*num_bodies), matvec = linear_operator_partial, dtype='float64')
      #Afull = np.empty((6*num_bodies, 6*num_bodies))
      #Id = np.identity(6*num_bodies)
      #for i in range(6*num_bodies):
	#if(i % 100 == 0):
	  #print i
	#Afull[:,i] = A * Id[:,i]
	
	
      #linear_operator_partial = partial(Lub.Diag_Lub_Mobility_Mult_RM, eta=eta, a=a, cutoff=4.5)
      #B = spla.LinearOperator((6*num_bodies, 6*num_bodies), matvec = linear_operator_partial, dtype='float64')
      #Bfull = np.empty((6*num_bodies, 6*num_bodies))
      #Id = np.identity(6*num_bodies)
      #for i in range(6*num_bodies):
	#if(i % 100 == 0):
	  #print i
	#Bfull[:,i] = B * Id[:,i]
	
      #print np.linalg.norm(Afull-Bfull)
      #plt.imshow(Afull-Bfull)
      #plt.show()
      
      

      F = np.random.rand(Nblobs*6)
      MF = Lub.Wall_Mobility_Mult(F, eta, a)
      Fb = np.concatenate((0.0*F,-1.0*F))
      print np.shape(Fb)
      UW, res = Lub.Lubrucation_alt_solve(MF, eta, a, 4.5, print_residual=True)
      #UW = Lub.Lubrucation_PC(Fb, eta, a, 4.5,its=1000)
      #UW, res = Lub.Lubrucation_solve(Fb, eta, a, 4.5)
      print(res)
      save_name = '/home/bs162/LubricationCode/RigidMultiblobsWallLubrication/GMRES_tests/Test_Configs_2/PC_'+bfname+'.gmres.txt'
      np.savetxt(save_name, res, delimiter='  ', fmt='%.15f')

