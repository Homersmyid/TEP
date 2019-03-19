# -*- coding: utf-8 -*-

from __future__ import division
from pyomo.environ import *
from pyomo.opt import SolverFactory
import numpy as np
import ruizC as RC
import ruizSub as s
import ruizMast as m
import ruizStep0 as szero

STOP = 2;							#How many iterations to quit after
startlines = True					#If possible lines at start
k = 0								#Iterations
LB = float("-inf")					#Upper Bound
UB = float("inf")					#Lower Bound

############################
#Step Zero Master
############################
#If there are lines at the start then find the cost of them
#This provides an absolute lower bound of the optimization
#Else the Lower Bound is zero	
if (startlines):
	#create step zero		
	imast = szero.mod.create_instance(RC.DATA)
		
	#Set x_star in step zero
	for x in RC.START_X_STAR:
		imast.x_star[x[0], x[1]] = x[2]
	
	#solve step zero
	zresults = szero.opt.solve(imast)
	print("\n\n***MASTER ZERO***\n\n")
	zresults.write()
	LB = value(imast.Obj)
	
	
	
	
	for v in imast.component_objects(Var, active=True):
		print ("Variable",v)
		varob = getattr(imast, str(v))
		for index in varob:
			print ("   ",index, varob[index].value)
	input()
	
	
else:
	LB = 0;

############################
#Step Zero Subproblem
############################
#Create subproblem
isub = s.mod.create_instance(RC.DATA)

#Set x_star in subproblem
for xi in imast.x:
	isub.x_star[xi] = int(round(value(imast.x[xi])))
		
#solve subproblem
sresults = s.opt.solve(isub)
print("\n\n***SUB ZERO***\n\n")
sresults.write()
UB = value(isub.Obj)



for v in isub.component_objects(Var, active=True):
	print ("Variable",v)
	varob = getattr(isub, str(v))
	for index in varob:
		print ("   ",index, varob[index].value)
#isub.pprint()
input()


############################
#Main Loop
############################
for k in range(1,STOP+1):
	
	#If UB and LB close enough, quit loop
	if abs(UB - LB) / UB <= RC.EPSILON:
		break
	
	########################
	#STEP K Master Problem
	######################## 
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
		imast.x_star[x[0], x[1]] = x[2]
	
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
	
	
	########################
	#STEP K Sub roblem
	########################
	#Create subproblem
	isub = s.mod.create_instance(RC.DATA)

	#Set x_star in sub
	for xi in imast.x:
		isub.x_star[xi] = int(round(value(imast.x[xi])))

	#solve subproblem
	sresults = s.opt.solve(isub)
	print('\n\nk:', k)
	print("*SUB***\n\n")
	sresults.write()
	UB = value(isub.Obj)
	
	for v in isub.component_objects(Var, active=True):
		print ("Variable",v)
		varob = getattr(isub, str(v))
		for index in varob:
			print ("   ",index, varob[index].value)
	input()
	

