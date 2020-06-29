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
from math import *
import numpy as np
import picos as sdp
import cvxopt as cvx

# global constants
display_float = '%.5f'
cvxopt = {'solver': 'cvxopt', 'solve_via_dual': False}


class BilinearTree:

    def __init__(self, bndVlow, bndVup, bndWlow, bndWup, treepos):
        self.bndVlow = np.array(bndVlow)
        self.bndVup = np.array(bndVup)
        self.bndWlow = np.array(bndWlow)
        self.bndWup = np.array(bndWup)
        self.children = []
        self.lower = None # lower bound - SDP value
        self.upper = None # upper bound - objective value
        self.solV = None
        self.solW = None
        self.treepos = treepos # string identifier


    def hyperarea(self):
        # calculate total area of the rectangles
        rk = len(self.bndVlow)
        area = 0.0
        for k in range(rk):
            area += (self.bndVup[k] - self.bndVlow[k]) * (self.bndWup[k] - self.bndWlow[k])
        return area


    def shrink_hyperrectangle(self, originalSDP, verb):
        # try to shrink hypercube
        bndSDP = originalSDP.copy()
        rk = len(self.bndVlow)

        oldvol = self.hyperarea()
        if oldvol == 0.0: return 0.0

        varV = bndSDP.get_variable('v')
        varW = bndSDP.get_variable('w')
        bndSDP.add_constraint(varV[0:rk] >= cvx.matrix(self.bndVlow))
        bndSDP.add_constraint(varV[0:rk] <= cvx.matrix(self.bndVup))
        bndSDP.add_constraint(varW[0:rk] >= cvx.matrix(self.bndWlow))
        bndSDP.add_constraint(varW[0:rk] <= cvx.matrix(self.bndWup))

        for k in range(rk):
            if (self.bndVup[k] - self.bndVlow[k]) * (self.bndWup[k] - self.bndWlow[k]) > 1.0e-6:
                try:
                    # X lower bound
                    bndSDP.set_objective('min', varV[k])
                    bndSDP.solve(verbose = verb-1, **cvxopt)
                    self.bndVlow[k] = bndSDP.obj_value()
                    # X upper bound
                    bndSDP.set_objective('max', varV[k])
                    bndSDP.solve(verbose = verb-1, **cvxopt)
                    self.bndVup[k] = bndSDP.obj_value()
                    # Y lower bound
                    bndSDP.set_objective('min', varW[k])
                    bndSDP.solve(verbose = verb-1, **cvxopt)
                    self.bndWlow[k] = bndSDP.obj_value()
                    # Y upper bound
                    bndSDP.set_objective('max', varW[k])
                    bndSDP.solve(verbose = verb-1, **cvxopt)
                    self.bndWup[k] = bndSDP.obj_value()

                except Exception as e:
                    print('sdp solver failed / maybe no solution')
                    continue

        return self.hyperarea() / oldvol


    def solve(self, baseSDP, a, b, Delta, verb):

        thisSDP = baseSDP.copy()
        rk = len(self.bndVlow)

        varV = thisSDP.get_variable('v')
        thisSDP.add_constraint(varV[0:rk] >= cvx.matrix(self.bndVlow))
        thisSDP.add_constraint(varV[0:rk] <= cvx.matrix(self.bndVup))

        varW = thisSDP.get_variable('w')
        thisSDP.add_constraint(varW[0:rk] >= cvx.matrix(self.bndWlow))
        thisSDP.add_constraint(varW[0:rk] <= cvx.matrix(self.bndWup))

        varR = thisSDP.get_variable('r')
        # construct SDP constraints specific for this branch
        thisSDP.add_constraint( (varV[0:rk]^cvx.matrix(self.bndWlow*Delta)) \
                + (varW[0:rk]^cvx.matrix(self.bndVlow*Delta)) \
                - cvx.matrix(self.bndVlow*self.bndWlow*Delta) \
                <= varR )
        thisSDP.add_constraint( (varV[0:rk]^cvx.matrix(self.bndWup*Delta)) \
                + (varW[0:rk]^cvx.matrix(self.bndVup*Delta)) \
                - cvx.matrix(self.bndVup*self.bndWup*Delta) \
                <= varR )

        # solve it
        self.lower = np.inf
        self.upper = np.inf

        try:
            thisSDP.solve(verbose = verb-1, **cvxopt)
        except Exception as e:
            print('sdp solver failed / maybe no solution')
            return False

        # extract solutions
        self.solV = np.array(thisSDP.get_valued_variable('v'))[:, 0]
        self.solW = np.array(thisSDP.get_valued_variable('w'))[:, 0]
        self.lower = thisSDP.obj_value()
        self.upper = np.sum(self.solV*a) + \
                np.sum(self.solW*b) + \
                np.sum(self.solV[0:rk]*Delta*self.solW[0:rk])

        return True


    def get_bounds(self, worstlower = np.inf, tobranch = None, bestupper = np.inf, best = None):

        if self.upper < bestupper:
            bestupper = self.upper
            best = self

        # we have reached a leaf
        if self.children == []:
            if self.lower < worstlower:
                worstlower = self.lower
                tobranch = self

        # we need to go further out the branches
        for child in self.children:
            (worstlower, tobranch, bestupper, best) = \
                child.get_bounds(worstlower, tobranch, bestupper, best)

        return (worstlower, tobranch, bestupper, best)


    def calculate_branch_index(self, rank, Delta, arot, brot):

        maxgap = -np.inf
        # only rank relevant here - for everything else the gap vanishes
        for h in range(0, rank):
            # shorthand notation
            l = self.bndVlow[h]; L = self.bndVup[h]
            m = self.bndWlow[h]; M = self.bndWup[h]
            v = self.solV[h]
            w = self.solW[h]
            gap = Delta[h]*(v*w - max(m*v + l*w - m*l, M*v + L*w - M*L))

            if gap > maxgap: idx = h; maxgap = gap

        return (idx, maxgap)


    def branch_here(self, allowedgap, rank, Delta, arot, brot, verb):

        (idx, gap) = self.calculate_branch_index(rank, Delta, arot, brot)
        # split the bounding hyperrectangle at the solution point
        newbndVlow = np.array(self.bndVlow)
        newbndVlow[idx] = self.solV[idx]
        newbndWlow = np.array(self.bndWlow)
        newbndWlow[idx] = self.solW[idx]
        newbndVup = np.array(self.bndVup)
        newbndVup[idx] = self.solV[idx]
        newbndWup = np.array(self.bndWup)
        newbndWup[idx] = self.solW[idx]

        if verb >= 2:
            print( ('branch at %s[%u]: ' + display_float \
                  + ' < ' + display_float + ' < ' + display_float + '; ' \
                  + display_float + ' < ' + display_float + ' < ' + display_float) \
                  % (self.treepos, idx, self.bndVlow[idx], self.solV[idx], self.bndVup[idx], \
                  self.bndWlow[idx], self.solW[idx], self.bndWup[idx]) )

        # descriptive string indicating position in tree
        path = '%s[%u]' % (self.treepos, idx)
        # create children going clockwise - starting at NE
        self.children.append(BilinearTree(newbndVlow, self.bndVup, newbndWlow, self.bndWup, path + '-ne'))
        self.children.append(BilinearTree(newbndVlow, self.bndVup, self.bndWlow, newbndWup, path + '-se'))
        self.children.append(BilinearTree(self.bndVlow, newbndVup, self.bndWlow, newbndWup, path + '-sw'))
        self.children.append(BilinearTree(self.bndVlow, newbndVup, newbndWlow, self.bndWup, path + '-nw'))
