# -*- coding: utf-8 -*-

from __future__ import division
from pyomo.environ import *
from pyomo.opt import SolverFactory
import ruizSub as s
import ruizMast as m
import ruizStep0 as szero

#Constants
SOLVER = "cplex"				#Solver to use
EPSILON = 10e-7					#Distance between up and lower bound to halt
DATA = "TEPTest.dat"			#Problem data file 
TOL = 10e-11					#Small Tolerance

#Lines existing at start
#In format ("Start Node", "End Node", "If connect (Bool)")


#No Lines At Start
START_X_STAR = [(1,2,0), (1,3,0), (1,4,0), (1,5,0), (1,6,0),
				(2,3,0), (2,4,0), (2,5,0), (2,6,0), (3,4,0),
				(3,5,0), (3,6,0), (4,5,0), (4,6,0), (5,6,0)]

STOP = 2;							#How many iterations to quit after
startlines = True					#If possible lines at start
k = 0								#Iterations
LB = float("-inf")					#Upper Bound
UB = float("inf")					#Lower Bound



############################
#Step Zero Subproblem
############################
#Create subproblem
isub = s.mod.create_instance(DATA)

#Set x_star in step zero
for x in START_X_STAR:
	isub.x_star[x[0], x[1]] = x[2]
		
#solve subproblem
sresults = s.opt.solve(isub)


print("\n\n***SUB ZERO***\n\n")

sresults.write()
'''
for v in isub.component_objects(Var, active=True):
	print ("Variable",v)
	varob = getattr(isub, str(v))
	for index in varob:
		print ("   ",index, varob[index].value)
'''
#isub.pprint()
input()









############################
#Main Loop
############################
for k in range(1,STOP+1):
	
	
	########################
	#STEP K Master Problem
	######################## 
	#create master problem
	if k == 1:			
		imast = m.mod.create_instance(DATA)

	m.mast_func(imast, isub.dem, isub.genpos, k)

	#solve master problem
	mresults = s.opt.solve(imast)
	print('\n\nk:', k)
	print("*MASTER*\n\n")
	mresults.write()
	LB = value(imast.Obj)

	
	for v in imast.component_objects(Var, active=True):
		print ("Variable",v)
		varob = getattr(imast, str(v))
		for index in varob:
			print ("   ",index, varob[index].value)
	
	#imast.pprint();
	input()

