# -*- coding: utf-8 -*-

#Constants
SOLVER = "cplex"				#Solver to use
EPSILON = 10e-7					#Distance between up and lower bound to halt
DATA = "data.dat"				#Problem data file 
TOL = 10e-11					#Small Tolerance

#Lines existing at start
#In format ("Start Node", "End Node", "If connect (Bool)")

'''
#No Lines At Start
START_X_STAR = [(1,2,0), (1,3,0), (1,4,0), (1,5,0), (1,6,0),
				(2,3,0), (2,4,0), (2,5,0), (2,6,0), (3,4,0),
				(3,5,0), (3,6,0), (4,5,0), (4,6,0), (5,6,0)]


				
'''
#Lines already existing at start: 	(1,2)&(1,4)&(1/5)&(2,3)&(2,4)&(3,5)
START_X_STAR = [(1,2,1), (1,3,0), (1,4,1), (1,5,1), (1,6,0),
				(2,3,1), (2,4,1), (2,5,0), (2,6,0), (3,4,0),
				(3,5,1), (3,6,0), (4,5,0), (4,6,0), (5,6,0)]

'''

#All Nodes Connected at Start
START_X_STAR = [(1,2,1), (1,3,1), (1,4,1), (1,5,1), (1,6,1),
				(2,3,1), (2,4,1), (2,5,1), (2,6,1), (3,4,1),
				(3,5,1), (3,6,1), (4,5,1), (4,6,1), (5,6,1)]
'''
