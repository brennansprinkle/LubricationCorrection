import argparse
import numpy as np
import scipy.linalg
import scipy.sparse.linalg as spla
import subprocess
import cPickle
from functools import partial
import sys
import time

# Add path to HydroGrid and import module
sys.path.append('../../HydroGrid/src/')

# Find project functions
found_functions = False
path_to_append = ''
while found_functions is False:
    try:
        import multi_bodies_functions
        from mobility import mobility as mb
        from quaternion_integrator.quaternion import Quaternion
        from quaternion_integrator.quaternion_integrator_multi_bodies import QuaternionIntegrator
        from quaternion_integrator.quaternion_integrator_rollers import QuaternionIntegratorRollers
        from body import body
        from read_input import read_input
        from read_input import read_vertex_file
        from read_input import read_clones_file
        from read_input import read_slip_file
        import utils

        try:
            import calculateConcentration as cc

            found_HydroGrid = True
        except ImportError:
            found_HydroGrid = False
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


def calc_slip(bodies, Nblobs):
    '''
    Function to calculate the slip in all the blobs.
    '''
    slip = np.empty((Nblobs, 3))
    offset = 0
    for b in bodies:
        slip_b = b.calc_slip()
        slip[offset:offset + b.Nblobs] = slip_b
        offset += b.Nblobs
    return slip


def get_blobs_r_vectors(bodies, Nblobs):
    '''
    Return coordinates of all the blobs with shape (Nblobs, 3).
    '''
    r_vectors = np.empty((Nblobs, 3))
    offset = 0
    for b in bodies:
        num_blobs = b.Nblobs
        r_vectors[offset:(offset + num_blobs)] = b.get_r_vectors()
        offset += num_blobs
    return r_vectors


def set_mobility_blobs(implementation):
    '''
    Set the function to compute the dense mobility
    at the blob level to the right implementation.
    The implementation in C++ is much faster than
    the one python; to use it the user should compile
    the file mobility/mobility_ext.cc.

    These functions return an array with shape
    (3*Nblobs, 3*Nblobs).
    '''
    # Implementations without wall
    if implementation == 'python_no_wall':
        return mb.rotne_prager_tensor
    elif implementation == 'C++_no_wall':
        return mb.boosted_infinite_fluid_mobility
    # Implementations with wall
    elif implementation == 'python':
        return mb.single_wall_fluid_mobility
    elif implementation == 'C++':
        return mb.boosted_single_wall_fluid_mobility


def set_mobility_vector_prod(implementation):
    '''
    Set the function to compute the matrix-vector
    product (M*F) with the mobility defined at the blob
    level to the right implementation.

    The implementation in pycuda is much faster than the
    one in C++, which is much faster than the one python;
    To use the pycuda implementation is necessary to have
    installed pycuda and a GPU with CUDA capabilities. To
    use the C++ implementation the user has to compile
    the file mobility/mobility_ext.cc.
    '''
    # Implementations without wall
    if implementation == 'python_no_wall':
        return mb.no_wall_fluid_mobility_product
    elif implementation == 'C++_no_wall':
        return mb.boosted_no_wall_mobility_vector_product
    elif implementation == 'pycuda_no_wall':
        return mb.no_wall_mobility_trans_times_force_pycuda
    # Implementations with wall
    elif implementation == 'python':
        return mb.single_wall_fluid_mobility_product
    elif implementation == 'C++':
        return mb.boosted_mobility_vector_product
    elif implementation == 'pycuda':
        return mb.single_wall_mobility_trans_times_force_pycuda


def calc_K_matrix(bodies, Nblobs):
    '''
    Calculate the geometric block-diagonal matrix K.
    Shape (3*Nblobs, 6*Nbodies).
    '''
    K = np.zeros((3 * Nblobs, 6 * len(bodies)))
    offset = 0
    for k, b in enumerate(bodies):
        K_body = b.calc_K_matrix()
        K[3 * offset:3 * (offset + b.Nblobs), 6 * k:6 * k + 6] = K_body
        offset += b.Nblobs
    return K


def calc_K_matrix_bodies(bodies, Nblobs):
    '''
    Calculate the geometric matrix K for
    each body. List of shape (3*Nblobs, 6*Nbodies).
    '''
    K = []
    for k, b in enumerate(bodies):
        K_body = b.calc_K_matrix()
        K.append(K_body)
    return K


def K_matrix_vector_prod(bodies, vector, Nblobs, K_bodies=None):
    '''
    Compute the matrix vector product K*vector where
    K is the geometrix matrix that transport the information from the
    level of describtion of the body to the level of describtion of the blobs.
    '''
    # Prepare variables
    result = np.empty((Nblobs, 3))
    v = np.reshape(vector, (len(bodies) * 6))

    # Loop over bodies
    offset = 0
    for k, b in enumerate(bodies):
        if K_bodies is None:
            K = b.calc_K_matrix()
        else:
            K = K_bodies[k]
        result[offset: offset + b.Nblobs] = np.reshape(np.dot(K, v[6 * k: 6 * (k + 1)]), (b.Nblobs, 3))
        offset += b.Nblobs
    return result


def K_matrix_T_vector_prod(bodies, vector, Nblobs, K_bodies=None):
    '''
    Compute the matrix vector product K^T*vector where
    K is the geometrix matrix that transport the information from the
    level of describtion of the body to the level of describtion of the blobs.
    '''
    # Prepare variables
    result = np.empty((len(bodies), 6))
    v = np.reshape(vector, (Nblobs * 3))

    # Loop over bodies
    offset = 0
    for k, b in enumerate(bodies):
        if K_bodies is None:
            K = b.calc_K_matrix()
        else:
            K = K_bodies[k]
        result[k: k + 1] = np.dot(K.T, v[3 * offset: 3 * (offset + b.Nblobs)])
        offset += b.Nblobs

    result = np.reshape(result, (2 * len(bodies), 3))
    return result


def linear_operator_rigid(vector, bodies, r_vectors, eta, a, K_bodies=None, *args, **kwargs):
    '''
    Return the action of the linear operator of the rigid body on vector v.
    The linear operator is
    |  M   -K|
    | -K^T  0|
    '''
    # Reserve memory for the solution and create some variables
    L = kwargs.get('periodic_length')
    Ncomp_blobs = r_vectors.size
    Nblobs = r_vectors.size / 3
    Ncomp_bodies = 6 * len(bodies)
    res = np.empty((Ncomp_blobs + Ncomp_bodies))
    v = np.reshape(vector, (vector.size / 3, 3))

    # Compute the "slip" part
    res[0:Ncomp_blobs] = mobility_vector_prod(r_vectors, vector[0:Ncomp_blobs], eta, a, *args, **kwargs)
    K_times_U = K_matrix_vector_prod(bodies, v[Nblobs: Nblobs + 2 * len(bodies)], Nblobs, K_bodies=K_bodies)
    res[0:Ncomp_blobs] -= np.reshape(K_times_U, (3 * Nblobs))

    # Compute the "-force_torque" part
    K_T_times_lambda = K_matrix_T_vector_prod(bodies, vector[0:Ncomp_blobs], Nblobs, K_bodies=K_bodies)
    res[Ncomp_blobs: Ncomp_blobs + Ncomp_bodies] = -np.reshape(K_T_times_lambda, (Ncomp_bodies))
    return res


@utils.static_var('mobility_bodies', [])
@utils.static_var('K_bodies', [])
@utils.static_var('M_factorization_blobs', [])
@utils.static_var('M_factorization_blobs_inv', [])
@utils.static_var('mobility_inv_blobs', [])
def build_block_diagonal_preconditioners_det_stoch(bodies, r_vectors, Nblobs, eta, a, *args, **kwargs):
    '''
    Build the deterministic and stochastic block diagonal preconditioners for rigid bodies.
    It solves exactly the mobility problem for each body
    independently, i.e., no interation between bodies is taken
    into account.

    If the mobility of a body at the blob
    level is M=L^T * L with L the Cholesky factor  we form the stochastic preconditioners

    P = inv(L)
    P_inv = L

    and the deterministic preconditioner
    N = (K.T * M^{-1} * K)^{-1}

    and return the functions to compute matrix vector products
    y = (P.T * M * P) * x
    y = P_inv * x
    y = N*F - N*K.T*M^{-1}*slip
    '''
    mobility_bodies = []
    K_bodies = []
    M_factorization_blobs = []
    M_factorization_blobs_inv = []
    mobility_inv_blobs = []

    if (kwargs.get('step') % kwargs.get('update_PC') == 0) or len(
            build_block_diagonal_preconditioners_det_stoch.mobility_bodies) == 0:
        print('making PC')
        # Loop over bodies
        for b in bodies:
            # 1. Compute blobs mobility
            M = b.calc_mobility_blobs(eta, a)
            # 2. Compute Cholesy factorization, M = L^T * L
            L, lower = scipy.linalg.cho_factor(M)
            L = np.triu(L)
            M_factorization_blobs.append(L.T)
            # 3. Compute inverse of L
            M_factorization_blobs_inv.append(scipy.linalg.solve_triangular(L, np.eye(b.Nblobs * 3), check_finite=False))
            # 4. Compute inverse mobility blobs
            mobility_inv_blobs.append(scipy.linalg.solve_triangular(L, scipy.linalg.solve_triangular(L, np.eye(
                b.Nblobs * 3), trans='T', check_finite=False), check_finite=False))
            # 5. Compute geometric matrix K
            K = b.calc_K_matrix()
            K_bodies.append(K)
            # 6. Compute body mobility
            mobility_bodies.append(
                np.linalg.pinv(np.dot(K.T, scipy.linalg.cho_solve((L, lower), K, check_finite=False))))

        # Save variables to use in next steps if PC is not updated
        build_block_diagonal_preconditioners_det_stoch.mobility_bodies = mobility_bodies
        build_block_diagonal_preconditioners_det_stoch.K_bodies = K_bodies
        build_block_diagonal_preconditioners_det_stoch.M_factorization_blobs = M_factorization_blobs
        build_block_diagonal_preconditioners_det_stoch.M_factorization_blobs_inv = M_factorization_blobs_inv
        build_block_diagonal_preconditioners_det_stoch.mobility_inv_blobs = mobility_inv_blobs
    else:
        # Use old values
        mobility_bodies = build_block_diagonal_preconditioners_det_stoch.mobility_bodies
        K_bodies = build_block_diagonal_preconditioners_det_stoch.K_bodies
        M_factorization_blobs = build_block_diagonal_preconditioners_det_stoch.M_factorization_blobs
        M_factorization_blobs_inv = build_block_diagonal_preconditioners_det_stoch.M_factorization_blobs_inv
        mobility_inv_blobs = build_block_diagonal_preconditioners_det_stoch.mobility_inv_blobs

    def block_diagonal_preconditioner(vector, bodies=None, mobility_bodies=None, mobility_inv_blobs=None, K_bodies=None,
                                      Nblobs=None, *args, **kwargs):
        '''
        Apply the block diagonal preconditioner.
        '''
        result = np.empty(vector.shape)
        offset = 0
        for k, b in enumerate(bodies):
            # 1. Solve M*Lambda_tilde = slip
            slip = vector[3 * offset: 3 * (offset + b.Nblobs)]
            Lambda_tilde = np.dot(mobility_inv_blobs[k], slip)
            # 2. Compute rigid body velocity
            F = vector[3 * Nblobs + 6 * k: 3 * Nblobs + 6 * (k + 1)]
            Y = np.dot(mobility_bodies[k], -F - np.dot(K_bodies[k].T, Lambda_tilde))
            # 3. Solve M*Lambda = (slip + K*Y)
            result[3 * offset: 3 * (offset + b.Nblobs)] = np.dot(mobility_inv_blobs[k], slip + np.dot(K_bodies[k], Y))
            # 4. Set result
            result[3 * Nblobs + 6 * k: 3 * Nblobs + 6 * (k + 1)] = Y
            offset += b.Nblobs
        return result

    block_diagonal_preconditioner_partial = partial(block_diagonal_preconditioner,
                                                    bodies=bodies,
                                                    mobility_bodies=mobility_bodies,
                                                    mobility_inv_blobs=mobility_inv_blobs,
                                                    K_bodies=K_bodies,
                                                    Nblobs=Nblobs)

    # Define preconditioned mobility matrix product
    def mobility_pc(w, bodies=None, P=None, r_vectors=None, eta=None, a=None, *args, **kwargs):
        result = np.empty_like(w)
        # Apply P
        offset = 0
        for k, b in enumerate(bodies):
            result[3 * offset: 3 * (offset + b.Nblobs)] = np.dot(P[k], w[3 * offset: 3 * (offset + b.Nblobs)])
            offset += b.Nblobs
        # Multiply by M
        result_2 = mobility_vector_prod(r_vectors, result, eta, a, *args, **kwargs)
        # Apply P.T
        offset = 0
        for k, b in enumerate(bodies):
            result[3 * offset: 3 * (offset + b.Nblobs)] = np.dot(P[k].T, result_2[3 * offset: 3 * (offset + b.Nblobs)])
            offset += b.Nblobs
        return result

    mobility_pc_partial = partial(mobility_pc, bodies=bodies, P=M_factorization_blobs_inv, r_vectors=r_vectors, eta=eta,
                                  a=a, *args, **kwargs)

    # Define inverse preconditioner P_inv
    def P_inv_mult(w, bodies=None, P_inv=None):
        offset = 0
        for k, b in enumerate(bodies):
            w[3 * offset: 3 * (offset + b.Nblobs)] = np.dot(P_inv[k], w[3 * offset: 3 * (offset + b.Nblobs)])
            offset += b.Nblobs
        return w

    P_inv_mult_partial = partial(P_inv_mult, bodies=bodies, P_inv=M_factorization_blobs)

    # Return preconditioner functions
    return block_diagonal_preconditioner_partial, mobility_pc_partial, P_inv_mult_partial


@utils.static_var('mobility_bodies', [])
@utils.static_var('K_bodies', [])
@utils.static_var('mobility_inv_blobs', [])
def build_block_diagonal_preconditioner(bodies, r_vectors, Nblobs, eta, a, *args, **kwargs):
    '''
    Build the block diagonal preconditioner for rigid bodies.
    It solves exactly the mobility problem for each body
    independently, i.e., no interation between bodies is taken
    into account.
    '''
    mobility_inv_blobs = []
    mobility_bodies = []
    K_bodies = []
    if (kwargs.get('step') % kwargs.get('update_PC') == 0) or len(
            build_block_diagonal_preconditioner.mobility_bodies) == 0:
        print('building PC')
        print(len(build_block_diagonal_preconditioner.mobility_bodies) == 0)
        # Loop over bodies
        for b in bodies:
            # 1. Compute blobs mobility and invert it
            M = b.calc_mobility_blobs(eta, a)
            # 2. Compute Cholesy factorization, M = L^T * L
            L, lower = scipy.linalg.cho_factor(M)
            L = np.triu(L)
            # 3. Compute inverse mobility blobs
            mobility_inv_blobs.append(scipy.linalg.solve_triangular(L, scipy.linalg.solve_triangular(L, np.eye(
                b.Nblobs * 3), trans='T', check_finite=False), check_finite=False))
            # 4. Compute geometric matrix K
            K = b.calc_K_matrix()
            K_bodies.append(K)
            # 5. Compute body mobility
            mobility_bodies.append(
                np.linalg.pinv(np.dot(K.T, scipy.linalg.cho_solve((L, lower), K, check_finite=False))))

        # Save variables to use in next steps if PC is not updated
        build_block_diagonal_preconditioner.mobility_bodies = mobility_bodies
        build_block_diagonal_preconditioner.K_bodies = K_bodies
        build_block_diagonal_preconditioner.mobility_inv_blobs = mobility_inv_blobs

    else:
        # Use old values
        mobility_bodies = build_block_diagonal_preconditioner.mobility_bodies
        K_bodies = build_block_diagonal_preconditioner.K_bodies
        mobility_inv_blobs = build_block_diagonal_preconditioner.mobility_inv_blobs

    def block_diagonal_preconditioner(vector, bodies=None, mobility_bodies=None, mobility_inv_blobs=None, K_bodies=None,
                                      Nblobs=None):
        '''
        Apply the block diagonal preconditioner.
        '''
        result = np.empty(vector.shape)
        offset = 0
        for k, b in enumerate(bodies):
            # 1. Solve M*Lambda_tilde = slip
            slip = vector[3 * offset: 3 * (offset + b.Nblobs)]
            Lambda_tilde = np.dot(mobility_inv_blobs[k], slip)

            # 2. Compute rigid body velocity
            F = vector[3 * Nblobs + 6 * k: 3 * Nblobs + 6 * (k + 1)]
            Y = np.dot(mobility_bodies[k], -F - np.dot(K_bodies[k].T, Lambda_tilde))

            # 3. Solve M*Lambda = (slip + K*Y)
            Lambda = np.dot(mobility_inv_blobs[k], slip + np.dot(K_bodies[k], Y))

            # 4. Set result
            result[3 * offset: 3 * (offset + b.Nblobs)] = Lambda
            result[3 * Nblobs + 6 * k: 3 * Nblobs + 6 * (k + 1)] = Y
            offset += b.Nblobs
        return result

    block_diagonal_preconditioner_partial = partial(block_diagonal_preconditioner,
                                                    bodies=bodies,
                                                    mobility_bodies=mobility_bodies,
                                                    mobility_inv_blobs=mobility_inv_blobs,
                                                    K_bodies=K_bodies,
                                                    Nblobs=Nblobs)
    return block_diagonal_preconditioner_partial


def block_diagonal_preconditioner(vector, bodies, mobility_bodies, mobility_inv_blobs, Nblobs):
    '''
    Block diagonal preconditioner for rigid bodies.
    It solves exactly the mobility problem for each body
    independently, i.e., no interation between bodies is taken
    into account.
    '''
    result = np.empty(vector.shape)
    offset = 0
    for k, b in enumerate(bodies):
        # 1. Solve M*Lambda_tilde = slip
        slip = vector[3 * offset: 3 * (offset + b.Nblobs)]
        Lambda_tilde = np.dot(mobility_inv_blobs[k], slip)

        # 2. Compute rigid body velocity
        F = vector[3 * Nblobs + 6 * k: 3 * Nblobs + 6 * (k + 1)]
        Y = np.dot(mobility_bodies[k], -F - np.dot(b.calc_K_matrix().T, Lambda_tilde))

        # 3. Solve M*Lambda = (slip + K*Y)
        Lambda = np.dot(mobility_inv_blobs[k], slip + np.dot(b.calc_K_matrix(), Y))

        # 4. Set result
        result[3 * offset: 3 * (offset + b.Nblobs)] = Lambda
        result[3 * Nblobs + 6 * k: 3 * Nblobs + 6 * (k + 1)] = Y
        offset += b.Nblobs
    return result


def build_stochastic_block_diagonal_preconditioner(bodies, r_vectors, eta, a, *args, **kwargs):
    '''
    Build block diagonal preconditioner to generate the noise
    for rigid bodies. If the mobility of a body at the blob
    level is M=L^T * L with L the Cholesky factor  we form the stochastic preconditioners

    P = inv(L)
    P_inv = L

    and return the functions to compute matrix vector products
    y = (P.T * M * P) * x
    y = P_inv * x
    '''
    P = []
    P_inv = []
    for b in bodies:
        # Compute blobs mobility for one body
        M = b.calc_mobility_blobs(eta, a)

        # 2. Compute Cholesy factorization, M = L^T * L
        L, lower = scipy.linalg.cho_factor(M)
        L = np.triu(L)
        P_inv.append(L.T)

        # Form preconditioners version P
        P.append(scipy.linalg.solve_triangular(L, np.eye(b.Nblobs * 3), check_finite=False))

    # Define preconditioned mobility matrix product
    def mobility_pc(w, bodies=None, P=None, r_vectors=None, eta=None, a=None, *args, **kwargs):
        result = np.empty_like(w)
        # Multiply by P.T
        offset = 0
        for k, b in enumerate(bodies):
            result[3 * offset: 3 * (offset + b.Nblobs)] = np.dot(P[k], w[3 * offset: 3 * (offset + b.Nblobs)])
            offset += b.Nblobs
        # Multiply by M
        result_2 = mobility_vector_prod(r_vectors, result, eta, a, *args, **kwargs)
        # Multiply by P
        offset = 0
        for k, b in enumerate(bodies):
            result[3 * offset: 3 * (offset + b.Nblobs)] = np.dot(P[k].T, result_2[3 * offset: 3 * (offset + b.Nblobs)])
            offset += b.Nblobs
        return result

    mobility_pc_partial = partial(mobility_pc, bodies=bodies, P=P, r_vectors=r_vectors, eta=eta, a=a, *args, **kwargs)

    # Define inverse preconditioner P_inv
    def P_inv_mult(w, bodies=None, P_inv=None):
        offset = 0
        for k, b in enumerate(bodies):
            w[3 * offset: 3 * (offset + b.Nblobs)] = np.dot(P_inv[k], w[3 * offset: 3 * (offset + b.Nblobs)])
            offset += b.Nblobs
        return w

    P_inv_mult_partial = partial(P_inv_mult, bodies=bodies, P_inv=P_inv)

    # Return preconditioner functions
    return mobility_pc_partial, P_inv_mult_partial




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
    n_steps = read.n_steps
    n_save = read.n_save
    n_relaxation = read.n_relaxation
    dt = read.dt
    eta = read.eta
    g = read.g
    a = read.blob_radius
    scheme = read.scheme
    output_name = read.output_name
    structures = read.structures
    structures_ID = read.structures_ID
    mobility_vector_prod = set_mobility_vector_prod(read.mobility_vector_prod_implementation)
    multi_bodies_functions.calc_blob_blob_forces = multi_bodies_functions.set_blob_blob_forces(
        read.blob_blob_force_implementation)
    multi_bodies_functions.calc_body_body_forces_torques = multi_bodies_functions.set_body_body_forces_torques(
        read.body_body_force_torque_implementation)

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
            b.mobility_blobs = set_mobility_blobs(read.mobility_blobs_implementation)
            b.ID = structures_ID[ID]
            # Calculate body length for the RFD
            if i == 0:
                b.calc_body_length()
            else:
                b.body_length = bodies[-1].body_length
            multi_bodies_functions.set_slip_by_ID(b, slip)
            # Append bodies to total bodies list
            bodies.append(b)
    bodies = np.array(bodies)

    # Set some more variables
    num_of_body_types = len(body_types)
    num_bodies = bodies.size
    Nblobs = sum([x.Nblobs for x in bodies])

    # Save bodies information
    with open(output_name + '.bodies_info', 'w') as f:
        f.write('num_of_body_types  ' + str(num_of_body_types) + '\n')
        f.write('body_names         ' + str(body_names) + '\n')
        f.write('body_types         ' + str(body_types) + '\n')
        f.write('num_bodies         ' + str(num_bodies) + '\n')
        f.write('num_blobs          ' + str(Nblobs) + '\n')

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
            b.mobility_blobs = set_mobility_blobs(read.mobility_blobs_implementation)
            b.ID = structures_ID[ID]
            # Calculate body length for the RFD
            if i == 0:
                b.calc_body_length()
            else:
                b.body_length = bodies[-1].body_length
            multi_bodies_functions.set_slip_by_ID(b, slip)
            # Append bodies to total bodies list
            bodies.append(b)
    bodies = np.array(bodies)

    bodies_base = []
    for k, b in enumerate(bodies):
        temp = np.zeros(3)
        temp[0] = bodies[k].location[0]
        temp[1] = bodies[k].location[1]
        temp[2] = bodies[k].location[2]
        bodies_base.append(temp)

    mob_coef = np.empty((2*200, 6))
    res_coef = np.empty((2*200, 6))
    evals = np.empty((200, 13))
    whole_matrices = np.empty((200, 12*12+1))
    d_count = -1
    for h in [199999]: #range(200)
        height = 1.0 + (h + 1.0) * 0.01 #
        dist = (h + 1.0) * 0.01
        x = height
        d_count += 1

        # Create rigid bodies

        for k, b in enumerate(bodies_base):
            bodies[k].location[0] = b[0]
            bodies[k].location[1] = b[1]
            bodies[k].location[2] = b[2]
            # print(bodies[k].location)

        #######################################################################
        for k, b in enumerate(bodies):
            bodies[k].location[2] -= 2.0
            factor = (2 + dist) / 3.0
            bodies[k].location *= factor
            bodies[k].location[2] += height
            print(bodies[k].location)
        #######################################################################
        r_vectors_blobs = get_blobs_r_vectors(bodies, Nblobs)
        slip = np.zeros(Nblobs*3)
        force_torque = np.zeros(num_bodies * 6)

        # Set right hand side
        System_size = Nblobs * 3 + num_bodies * 6
        RHS = np.reshape(np.concatenate([slip, -force_torque]), (System_size))

        # Set linear operators
        linear_operator_partial = partial(linear_operator_rigid, bodies=bodies, r_vectors=r_vectors_blobs, eta=read.eta,
                                          a=read.blob_radius)
        A = spla.LinearOperator((System_size, System_size), matvec=linear_operator_partial, dtype='float64')

        # Set preconditioner
        if (d_count % 20 == 0):
            mobility_inv_blobs = []
        mobility_bodies = np.empty((len(bodies), 6, 6))
        # Loop over bodies
        for k, b in enumerate(bodies):
            # 1. Compute blobs mobility and invert it
            if (d_count % 20 == 0):
                M = b.calc_mobility_blobs(read.eta, read.blob_radius)
                M_inv = np.linalg.inv(M)
                mobility_inv_blobs.append(M_inv)
            # 2. Compute body mobility
            N = b.calc_mobility_body(read.eta, read.blob_radius, M_inv=M_inv)
            mobility_bodies[k] = N

        # 4. Pack preconditioner
        PC_partial = partial(block_diagonal_preconditioner, bodies=bodies, mobility_bodies=mobility_bodies, \
                             mobility_inv_blobs=mobility_inv_blobs, Nblobs=Nblobs)
        PC = spla.LinearOperator((System_size, System_size), matvec=PC_partial, dtype='float64')

        Lub_Mob = np.zeros((12, 12))
        start = time.time()
        
        class gmres_counter(object):
	  def __init__(self, disp=True):
	      self._disp = disp
	      self.niter = 0
	  def __call__(self, rk=None):
	      self.niter += 1
	      if self._disp:
		  print('iter %3i\trk = %s' % (self.niter, str(rk)))
        
        
        for p in range(12):
            # print(p)
            force_torque *= 0.0
            force_torque[p] = 1.0
            # print force_torque
            # Set right hand side
            RHS = np.reshape(np.concatenate([slip, -force_torque]), (System_size))
            # Solve preconditioned linear system # callback=make_callback()
            counter = gmres_counter()
            (sol_precond, info_precond) = utils.gmres(A, RHS, tol=read.solver_tolerance, M=PC, maxiter=10000, restart=1000, callback=counter)
            # Extract velocities and constraint forces on blobs
            velocity = sol_precond[3 * Nblobs: 3 * Nblobs + 6 * num_bodies]
            lambda_blobs = sol_precond[0: 3 * Nblobs]

            Lub_Mob[:, p] = velocity

        end = time.time()
        print(end - start)
        P = np.kron(np.array([[1., 0., 0., 0.], [0., 0., 1., 0.], [0., 1., 0., 0.], [0., 0., 0., 1.]]), np.eye(3))
        #Mob_MB = np.dot(P, np.dot(Lub_Mob, P))
        #Lub_Resistance_MB = np.linalg.pinv(Lub_Mob)
        #Resistance_MB = np.linalg.pinv(Mob_MB)
        
        
        
        # print(Resistance_MB)
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
	whole_matrix = Lub_Mob.flatten()
	whole_matrix = np.insert(whole_matrix, 0, x, axis=0)
	whole_matrices[d_count, :] = whole_matrix
	#idx = eigenValues.argsort()[::-1]   
	#eigenValues = eigenValues[idx]
	#evals[d_count, :] = eigenValues

    #mob_name = '/home/bs162/LubricationCode/RigidMultiblobsWallLubrication/multi_bodies/examples/Tetra_Multiblobs/Pair_Mob_Data/FarFEO_0.25_mob_coefs_642.txt' #FarFW_ FarFEO_
    #np.savetxt(mob_name, mob_coef, delimiter='  ', fmt='%.15f')
    #res_name = '/home/bs162/LubricationCode/RigidMultiblobsWallLubrication/multi_bodies/examples/Tetra_Multiblobs/Pair_Mob_Data/FarFEO_0.25_res_coefs_642.txt'
    #np.savetxt(res_name, res_coef, delimiter='  ')
    
    print Lub_Mob
    #mob_name = '/home/bs162/LubricationCode/RigidMultiblobsWallLubrication/multi_bodies/Pair_Mob_Data/mob_matrix_2562_e_20.txt'
    #np.savetxt(mob_name, whole_matrices, delimiter='  ', fmt='%.15f')

