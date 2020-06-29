"""
© 2018 Stefan Huber, Robert König, Marco Tomamichel

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
from tree import *
import numpy as np
import picos as sdp
import cvxopt as cvx

# some global parameters
STD_MAXITER = 100000
STD_ALLOWEDGAP = 1e-5


class BilinearProblem:

    def __init__(self):
        self.dim = [0,0]
        self.rank = 0
        self.Delta = []
        self.S = [[]]
        self.T = [[]]
        self.Q = [[]]
        self.a = []
        self.b = []
        self.maximize = False
        self.SDP = None
        self.varV = None
        self.varW = None
        self.solX = None
        self.solY = None
        self.bndVlow = []
        self.bndVup = []
        self.bndWlow = []
        self.bndWup = []
        self.tree = None
        self.lower = 0
        self.upper = 0
        self.bestleaf = None
        self.tosolve = []
        self.iter = 0


    def init_vector_form(self, Q, a = [], b = [], maximize = False, verb = 1):
        assert(self.SDP == None) # only once

        # set internals
        self.Q = np.array(Q)
        dim = self.dim = self.Q.shape

        if a == []: self.a = np.zeros(dim[0])
        else: self.a = np.array(a)
        assert(len(self.a) == dim[0])

        if b == []: self.b = np.zeros(dim[1])
        else: self.b = np.array(b)
        assert(len(self.b) == dim[1])

        self.maximize = maximize
        self.SDP = baseSDP = sdp.Problem()

        # Q = S D T with S, T unitary, D diagonal (represented as vector)
        (self.S, Delta, self.T) = np.linalg.svd(Q, full_matrices=True)

        # cleanup and calculate rank
        norm = Delta[0]
        self.rank = min(dim)
        for idx in range(1, min(dim)):
            if Delta[idx] < 1e-10*norm:
                self.rank = idx
                break
        self.Delta = Delta[0:self.rank]

        if verb >= 1:
            print('initialize vector problem: dim = (%u, %u), rank = %u' % (dim[0], dim[1], self.rank))
            print('Delta =', self.Delta)
            print('S =', self.S)
            print('T =', self.T)

        # setup sdp variables, w = T y, v = S^* x
        self.varV = baseSDP.add_variable('v', dim[0])
        self.varW = baseSDP.add_variable('w', dim[1]) # helper variable

    def add_constraint(self, constr):
        self.SDP.add_constraint(constr)

    def varX(self):
        return cvx.matrix(self.S)*self.varV

    def varY(self):
        return cvx.matrix(self.T).H*self.varW


    def solve(self, verb = 1, maxiter = STD_MAXITER, allowedgap = STD_ALLOWEDGAP):
        assert(self.maximize == False); # currently not supported
        # check if this is the first time
        assert(self.tosolve == [])

        # step 1: determine hyperrectangle of interest (for x and w variables)
        if self.bndVlow == []:
            self.init_hyperrectangle(verb)
        # initialize tree
        self.tree = BilinearTree(self.bndVlow, self.bndVup, \
                                 self.bndWlow, self.bndWup, 'rt')

        self.tosolve.append(self.tree)
        self.iter = 0

        self.solvemore(verb, maxiter, allowedgap)


    def solvemore(self, verb = 1, maxiter = STD_MAXITER, allowedgap = STD_ALLOWEDGAP):

        # are we done already?
        if self.iter >= maxiter+1:
            print('no more iterations!')
            return

        # step 2: initialize the SDP
        baseSDP = self.SDP.copy()
        varR = baseSDP.add_variable('r', self.rank)

        # rotate a and b matrices
        arot = np.array(np.matrix(self.a)*np.matrix(self.S))[0,:]
        brot = np.array(np.matrix(self.b)*np.matrix(self.T).H)[0,:]
        baseSDP.set_objective('min', (1|varR) \
                              + (cvx.matrix(arot)|self.varV) \
                              + (cvx.matrix(brot)|self.varW) )

        # step 3: iterate solve and branch
        for i in range(self.iter, maxiter+1):

            for leaf in self.tosolve:
                leaf.solve(baseSDP, arot, brot, self.Delta, verb)

            # calculate bounds
            (lower, leaf, upper, best) = self.tree.get_bounds()
            if verb >= 1:
                print( ('iteration: %6u, best: ' + display_float + ', lower: ' + display_float + ', gap: ' + display_float) % \
                      (i, upper, lower, upper - lower) )

            # check if we are done
            if upper - lower < allowedgap:
                # we are done
                print('solution found! terminating...')
                break

            assert(leaf != None)
            leaf.branch_here(allowedgap, self.rank, self.Delta, arot, brot, verb)

            self.tosolve = leaf.children

        # next round of iterations
        self.iter = maxiter + 1
        # copy the result here
        self.worstleaf = leaf
        self.bestleaf = best
        self.solX = np.array(np.matrix(self.S)*np.matrix(best.solV).H)[:,0]
        self.solY = np.array(np.matrix(self.T).H*np.matrix(best.solW).H)[:,0]

        self.lower = lower
        self.upper = upper


    def init_hyperrectangle(self, verb = 1):
        assert(self.bndVlow == []); # only once

        info = 'bounding hyperrectangle: '

        bndSDP = self.SDP.copy()
        varV = bndSDP.get_variable('v')
        varW = bndSDP.get_variable('w')
        for k in range(self.rank):
            # X lower bound
            bndSDP.set_objective('min', varV[k])
            bndSDP.solve(verbose = verb-1, **cvxopt)
            self.bndVlow.append(bndSDP.obj_value())
            # X upper bound
            bndSDP.set_objective('max', varV[k])
            bndSDP.solve(verbose = verb-1, **cvxopt)
            self.bndVup.append(bndSDP.obj_value())
            # Y lower bound
            bndSDP.set_objective('min', varW[k])
            bndSDP.solve(verbose = verb-1, **cvxopt)
            self.bndWlow.append(bndSDP.obj_value())
            # Y upper bound
            bndSDP.set_objective('max', varW[k])
            bndSDP.solve(verbose = verb-1, **cvxopt)
            self.bndWup.append(bndSDP.obj_value())
            # dispaly results
            if verb >= 1:
                print( ('%s ' + display_float + ' < v[%u] < ' + display_float \
                        + ', ' + display_float + ' < w[%u] < ' + display_float) % \
                        (info, self.bndVlow[-1], k, self.bndVup[-1], self.bndWlow[-1], k, self.bndWup[-1]) )
                info = '                         '
