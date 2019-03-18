# -*- coding: utf-8 -*-

from __future__ import division
from pyomo.environ import *
from pyomo.opt import SolverFactory
import numpy as np
import ruizC as RC
import ruizSub as s
import ruizMast as m


#Set Up
k = 0								#Iterations
LB = float("-inf")					#Upper Bound
UB = float("inf")					#Lower Bound

#Create subproblem
isub = s.sub.create_instance(RC.DATA)

#Set x_star in subproblem
for x in RC.START_X_STAR:
	isub.x_star[x[0], x[1]] = x[2]

#solve subproblem
results = s.opt.solve(isub)
print("\n***FIRSTSUB***\n\n")
isub.pprint()
results.write()
input()
	
for k in range(1,2+1):
	print('\n\nk:')
	print(k)
	
	#create master problem			
	imast = m.mod.create_instance(RC.DATA)

	#set demand in master
	for d in isub.dem:
		imast.dem[d] = value(isub.dem[d])

	#set possible generation in master
	for g in isub.genpos:
		imast.genpos[g] = value(isub.genpos[g])
		
	#Set x_star in master
	for x in RC.START_X_STAR:
		isub.x_star[x[0], x[1]] = x[2]
	
	#solve master problem
	results = s.opt.solve(imast)
	print("\n\n***MASTER***\n\n")
	imast.pprint()
	results.write()
	input()
	
	#Create subproblem
	isub = s.sub.create_instance(RC.DATA)

	#Set x_star in sub
	for xi in imast.x:
		isub.x_star[xi] = int(round(value(imast.x[xi])))

	#solve subproblem
	results = s.opt.solve(isub)
	print("\n***SUB***\n\n")
	isub.pprint()
	results.write()
	input()

