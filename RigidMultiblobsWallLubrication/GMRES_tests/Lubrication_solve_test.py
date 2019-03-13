import argparse
import numpy as np
import scipy.linalg as sp
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
  
  

  #num_bodies = 2
  
  #F_1 = np.zeros(6*num_bodies)
  #F_2 = np.zeros(6*num_bodies)
  #eps = 0.01
  
  
  #Rx11a = []
  #Rx12a = []
  #Ry11a = []
  #Ry12a = []
  #Ry11b = []
  #Ry12b = []
  #Rx11c = []
  #Rx12c = []
  #Ry11c = []
  #Ry12c = []
  
  #P = np.kron(np.array([[1., 0., 0., 0.],[0., 0., 1., 0.],[0., 1., 0., 0.],[0., 0., 0., 1.]]), np.eye(3))
  #sp_count = 0
  #for z_pt in [1.1*a, 1.5*a, 2.49*a, 2.51*a, 3.0*a, 10000.0*a]: #
    #sp_count += 1
    #x = [];
    #min_e_val = []
    #cutoff = 4.5
    #for j in range(299):
      #x_pt = (5.0-j*eps)*a#3.05
      ##z_pt = 1.1*a#(3.0-j*eps)
      #new_bodies = []
      #new_bodies.append(bodies[0])
      #new_bodies.append(bodies[1])
      #new_bodies[0].location = np.array([0.0, 0.0, z_pt])
      #new_bodies[1].location = np.array([x_pt, 0.0, z_pt])
      #Lub.bodies = new_bodies
      #x.append(x_pt)
      #print(j)
      #Lub_Cor_Mob = np.zeros((12,12))
      #RP_Mob = np.zeros((12,12))
      #for i in range(12):
	#F_2 = 0.0*F_2
	#F_2[i] = 1.0
	#F = np.concatenate((F_1,F_2))
	#UW, res = Lub.Lubrucation_alt_solve(F_2, eta, a, cutoff)
	#UW_R = Lub.pair_resistance_mult_tree(F_2, eta, a, cutoff) - Lub.pair_resistance_mult_MB_tree(F_2, eta, a, cutoff)
	#Lub_Cor_Mob[:,i] = UW#[12::]
	#RP_Mob[:,i] = UW_R
      #Lub_Cor_Mob = np.dot(P,np.dot(Lub_Cor_Mob,P))
      #Delta_R = np.dot(P,np.dot(RP_Mob,P))
      #evals = sp.eigh(Delta_R, eigvals_only=True)
      ##evals = sp.eigh(Lub_Cor_Mob, eigvals_only=True)
      #min_e_val.append(min(evals))
      
      #Rx11a.append(Lub_Cor_Mob[0,0])
      #Rx12a.append(Lub_Cor_Mob[0,3])
      #Ry11a.append(Lub_Cor_Mob[1,1])
      #Ry12a.append(Lub_Cor_Mob[1,4])
      #Ry11b.append((-Lub_Cor_Mob[4,11]))
      #Ry12b.append((-Lub_Cor_Mob[1,11]))
      #Rx11c.append(Lub_Cor_Mob[9,9])
      #Rx12c.append(Lub_Cor_Mob[9,6])
      #Ry11c.append(Lub_Cor_Mob[11,11])
      #Ry12c.append(Lub_Cor_Mob[11,8])
      
      
    ##print min_e_val[0] 
    ##plt.rcParams.update({'font.size': 20})
    ##plt.plot(x,min_e_val,'-ro')
    ##plt.title('smallest eigenvalue for Delta R with 2 particlels and no wall, assymptote = ' + str(min_e_val[0]), fontsize=30)
    ##plt.ylabel('smallest eigenvalue', fontsize=30)
    ##plt.xlabel('center distance', fontsize=30)
    ##plt.show()
      
    #ax = plt.subplot(2,3,sp_count) 
    #ax.set_title("height = " + str(z_pt), fontsize=30)
    #if sp_count > 3:
      #plt.xlabel('center distance', fontsize=30)
    #plt.rcParams.update({'font.size': 20})
    #if sp_count == 1:
      #plt.ylabel('smallest eigenvalue', fontsize=30)
    #plt.plot(x,min_e_val,'-ro')
  #plt.show()
  
  
  
  
  
  #plt.subplot(5,2,1)
  #plt.plot(x,Rx11a,'-ro')
  #plt.subplot(5,2,2)
  #plt.plot(x,Rx12a,'-ro')
  #plt.subplot(5,2,3)
  #plt.plot(x,Ry11a,'-ro')
  #plt.subplot(5,2,4)
  #plt.plot(x,Ry12a,'-ro')
  #plt.subplot(5,2,5)
  #plt.plot(x,Ry11b,'-ro')
  #plt.subplot(5,2,6)
  #plt.plot(x,Ry12b,'-ro')
  #plt.subplot(5,2,7)
  #plt.plot(x,Rx11c,'-ro')
  #plt.subplot(5,2,8)
  #plt.plot(x,Rx12c,'-ro')
  #plt.subplot(5,2,9)
  #plt.plot(x,Ry11c,'-ro')
  #plt.subplot(5,2,10)
  #plt.plot(x,Ry12c,'-ro')
  #plt.show()
  
  F = np.zeros(6)
  eps = 0.01

  x = [];
  min_e_val = []
  cutoff = 4.5
  for j in range(150):
    z_pt = (2.5-j*eps)
    new_bodies = []
    new_bodies.append(bodies[0])
    new_bodies[0].location = np.array([0.0, 0.0, z_pt])
    Lub.bodies = new_bodies
    x.append(z_pt)
    print(z_pt)
    Wall_Mob = Lub.Lub_Resist_Wall(z_pt, eta, a) - Lub.Lub_Resist_Wall_MB(z_pt, eta, a)
    evals = sp.eigh(Wall_Mob, eigvals_only=True)
    min_e_val.append(min(evals))
    
    
  I = [ n for n,i in enumerate(min_e_val) if i>0.0 ][0]
  print min_e_val[0]
  print x[I]
  plt.rcParams.update({'font.size': 20})
  plt.plot(x,min_e_val,'-ro')
  plt.title(r"smallest eigenvalue for (6,6) wall correction, assymptote = " + str(min_e_val[0]) + "\n" + "first height with negative eigen value: " + str(x[I-1]), fontsize=30)
  plt.ylabel('smallest eigenvalue', fontsize=30)
  plt.xlabel('wall distance to particle center', fontsize=30)
  plt.show()
  

  
  