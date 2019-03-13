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
#import matplotlib.pyplot as plt

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
        print
        'Creating structures = ', structure[1]
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

    # Res_Wall_Coeffs = Lub.Single_Blob_Resistance(eta, a, 1502, 0.0009998)
    mob_coef = np.empty((2*500, 6))
    res_coef = np.empty((2*500, 6))
    evals = np.empty((500, 13))
    whole_matrices = np.empty((500, 12*12+1))
    d_count = -1
    for h in range(500): #
        height = 1.0 + (h + 1.0) * 0.01
        dist = (h + 1.0) * 0.01
        x = height
        d_count += 1
        new_bodies = []
        for k, b in enumerate(bodies):
            n_b = copy.deepcopy(b)
            new_bodies.append(n_b)
            new_bodies[k].location[2] -= 2.0
            factor = (2 * a + dist * a) / 3.0
            new_bodies[k].location *= factor
            new_bodies[k].location[2] += height
            print(new_bodies[k].location)
        Lub.bodies = new_bodies
        F = np.zeros(2 * 6)
        F_s = np.zeros(2 * 2 * 6)
        Lub_Mob = np.zeros((12, 12))
        start = time.time()
        for p in range(12):
            # print(p)
            F *= 0.0
            G = copy.deepcopy(F)
            F[p] = 1.0
            F_s *= 0.0
            F_s = np.concatenate((0.0*G,-1.0*F))
            #derp = Lub.pair_resistance_mult_MB_tree(F, eta, a, 4.5) 
            #Itil = Lub.Wall_Mobility_Mult(derp, eta, a)
            UW, res = Lub.Lubrucation_alt_solve(F, eta, a, 4.5) #Lub.Lubrucation_solve(F_s, eta, a, 4.5) #
            Lub_Mob[:, p] = UW#[2 * 6::] #Itil#
  
	#print Lub_Mob
        end = time.time()
        print(end - start)
        P = np.kron(np.array([[1., 0., 0., 0.], [0., 0., 1., 0.], [0., 1., 0., 0.], [0., 0., 0., 1.]]), np.eye(3))
        
        #TODO CHANGE ME BACK FOR FARFW
        Mob_MB = np.dot(P, np.dot(Lub_Mob, P))
        
        
        #print(Mob_MB)
        #print(Lub_Mob)
        Lub_Resistance_MB = np.linalg.pinv(Lub_Mob)
        Resistance_MB = np.linalg.pinv(Mob_MB)
        # print(Resistance_MB)
        
        ## HACK USE THIS WHEN ONLY DIST BETWEEN PARTICLES IS VARIED
        #res_coef_j_11 = np.array([x, Resistance_MB[0, 0], Resistance_MB[1, 1], Resistance_MB[4, 11], Resistance_MB[9, 9]
                                     #, Resistance_MB[11, 11]])
        #res_coef_j_12 = np.array([x, Resistance_MB[0, 3], Resistance_MB[1, 4], Resistance_MB[1, 11], Resistance_MB[9, 6]
                                     #, Resistance_MB[11, 8]])
        #res_coef[2 * d_count, :] = res_coef_j_11
        #res_coef[2 * d_count + 1, :] = res_coef_j_12

        #mob_coef_j_11 = np.array([x, Mob_MB[0, 0], Mob_MB[1, 1], Mob_MB[4, 11], Mob_MB[9, 9], Mob_MB[11, 11]])
        #mob_coef_j_12 = np.array([x, Mob_MB[0, 3], Mob_MB[1, 4], Mob_MB[1, 11], Mob_MB[9, 6], Mob_MB[11, 8]])
        #mob_coef[2 * d_count, :] = mob_coef_j_11
        #mob_coef[2 * d_count + 1, :] = mob_coef_j_12
        
        
        #res_coef_j_11 = np.array([x, Lub_Resistance_MB[2,2], Lub_Resistance_MB[0,0], Lub_Resistance_MB[4,0], Lub_Resistance_MB[5,5], Lub_Resistance_MB[3,3]])
        #res_coef_j_12 = np.array([x, Resistance_MB[0,3], Resistance_MB[1,4], Resistance_MB[1,11], Resistance_MB[9,6], Resistance_MB[11,8]])
        #res_coef[2 * d_count, :] = res_coef_j_11
        #res_coef[2 * d_count + 1, :] = res_coef_j_12

        #mob_coef_j_11 = np.array([x, Lub_Mob[2,2], Lub_Mob[0,0], Lub_Mob[4,0], Lub_Mob[5,5], Lub_Mob[3,3]])
        #mob_coef_j_12 = np.array([x, Mob_MB[0,3], Mob_MB[1,4], Mob_MB[1,11], Mob_MB[9,6], Mob_MB[11,8]])
        #mob_coef[2 * d_count, :] = mob_coef_j_11
        #mob_coef[2 * d_count + 1, :] = mob_coef_j_12
        
        #eigenValues, eigenVectors = scipy.linalg.eig(Lub_Mob)
	#eigenValues = np.insert(eigenValues, 0, x, axis=0)
	#idx = eigenValues.argsort()[::-1]   
	#eigenValues = eigenValues[idx]
	#evals[d_count, :] = eigenValues
	
	whole_matrix = Lub_Mob.flatten()
	whole_matrix = np.insert(whole_matrix, 0, x, axis=0)
	whole_matrices[d_count, :] = whole_matrix
	
    #mob_name = '/home/bs162/LubricationCode/RigidMultiblobsWallLubrication/multi_bodies/Pair_Mob_Data/FarFEO_0.25_mob_coefs.txt' #FarFEO_3_
    #np.savetxt(mob_name, mob_coef, delimiter='  ', fmt='%.15f')
    #res_name = '/home/bs162/LubricationCode/RigidMultiblobsWallLubrication/multi_bodies/Pair_Mob_Data/FarFEO_0.25_res_coefs.txt' #FarFEO_3_
    #np.savetxt(res_name, res_coef, delimiter='  ')
  
    mob_name = '/home/bs162/LubricationCode/RigidMultiblobsWallLubrication/multi_bodies/Pair_Mob_Data/mob_matrix_PC_test.txt'
    np.savetxt(mob_name, whole_matrices, delimiter='  ', fmt='%.15f')
