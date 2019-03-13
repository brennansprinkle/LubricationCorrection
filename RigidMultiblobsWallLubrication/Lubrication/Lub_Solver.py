'''
Class to handle Lubrication solve
'''
import numpy as np
import scipy.spatial as spatial
import scipy.sparse.linalg as spla
from scipy.linalg import cholesky
from scipy.linalg import solve_triangular
from functools import partial
import copy
import inspect
import time
import utils
import sys
from LLC import LLC

class Lub_Solver(object):
  '''
  Class to handle Lubrication solve
  '''  
  def __init__(self, bodies):
    '''
    Constructor. Take arguments like ...
    '''
    # Location as np.array.shape = 3
    self.bodies = bodies
    self.mob_scalars_MB = None
    self.res_scalars_wall_MB_1 = None
    self.mob_scalars_wall_MB = None
    self.mob_scalars_WS = None
    self.mob_scalars_JO = None
    self.MB_Fn = None
    self.MB_wall_Fn = None
    self.MB_res_wall_Fn = None
    self.WS_Fn = None
    
    # for PC
    self.Lub_Chol_Blocks = None
    self.Mob_Inv_Blocks = None
    
    #if domain is a 'single_wall'
    self.mobility_trans_times_force_wall = mob.single_wall_mobility_trans_times_force_numba
    self.mobility_trans_times_torque_wall = mob.single_wall_mobility_trans_times_torque_numba
    self.mobility_rot_times_force_wall = mob.single_wall_mobility_rot_times_force_numba
    self.mobility_rot_times_torque_wall = mob.single_wall_mobility_rot_times_torque_numba
    #if domain is a 'no_wall'
    self.mobility_trans_times_force = mob.no_wall_mobility_trans_times_force_numba
    self.mobility_trans_times_torque = mob.no_wall_mobility_trans_times_torque_numba
    self.mobility_rot_times_force = mob.no_wall_mobility_rot_times_force_numba
    self.mobility_rot_times_torque = mob.no_wall_mobility_rot_times_torque_numba
  
  def load_WS_coefficient_interp_data(self):
    self.mob_scalars_WS = np.loadtxt("./Resistance_Coefs/mob_scalars_WS.txt") #res_scalars_WS.txt