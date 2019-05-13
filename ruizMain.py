# -*- coding: utf-8 -*-

from __future__ import division
from pyomo.environ import *
from pyomo.opt import SolverFactory
from pyomo.opt import SolverStatus, TerminationCondition
import ruizC as RC
import ruizSub as s
import ruizMast as m
import time
start = time.time()

STOP = 8							#How many iterations to quit after
startlines = True					#If possible lines at start
LB = float("-inf")					#Upper Bound
UB = float("inf")					#Lower Bound

#Lines At Start
START_X_STAR = [(1,2,1), (1,3,1), (1,4,0), (1,5,0), (1,6,0),
				(2,3,0), (2,4,0), (2,5,0), (2,6,0), (3,4,0),
				(3,5,0), (3,6,0), (4,5,1), (4,6,0), (5,6,0)]

############################
#Step Zero Master
############################
#If there are lines at the start then find the cost of them
#This provides an absolute lower bound of the optimization
#Else the Lower Bound is zero	
if (startlines):
	#create step zero		
	imast = m.mod.create_instance(RC.DATA)
		
	#Set x_star in step zero
	for x in START_X_STAR:
		imast.x_star[x[0], x[1]] = x[2]
	
	#solve step zero
	zresults = m.opt.solve(imast)
	LB = value(imast.Obj)

	
	print("\n\n***MASTER ZERO***\n\n")
	zresults.write()
	for v in imast.component_objects(Var, active=True):
		print ("Variable",v)
		varob = getattr(imast, str(v))
		for index in varob:
			print ("   ",index, varob[index].value)
	
	input()
	
	
else:
	LB = 0

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
UB = value(isub.Obj)


print("\n\n***SUB ZERO***\n\n")
sresults.write()
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
	if (UB - LB) / UB <= RC.EPSILON:
		break
	
	########################
	#STEP K Master Problem
	######################## 
	m.mast_func(imast, isub.dem, isub.genpos, START_X_STAR, k)

	#solve master problem
	mresults = m.opt.solve(imast)
	LB = value(imast.Obj)
	
	print('\n\nk:', k)
	print("*MASTER*\n\n")
		
	mresults.write()	
	for v in imast.component_objects(Var, active=True):
		print ("Variable",v)
		varob = getattr(imast, str(v))
		for index in varob:
			print ("   ",index, varob[index].value)
	#imast.pprint()	
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
	for v in isub.component_objects(Var, active=True):
		print ("Variable",v)
		varob = getattr(isub, str(v))
		for index in varob:
			print ("   ",index, varob[index].value)
	#isub.pprint()
	input()
	
	
	#store new upper bound if it is < previous upper bound
	if value (isub.Obj) <= UB:
		print("UPDATE UB")
		UB = value(isub.Obj)
	
	print("XXX")
	print(UB)
	print(LB)
	print(UB - LB)
	print((UB-LB)/ UB)
	print("XXX")
	
	end = time.time()
	print("Time:")
	print(end - start)
	
	
'''
results = opt.solve(instance) # Solving a model instance  
instance.load(results) # Loading solution into results object

if (results.solver.status == SolverStatus.ok) and (results.solver.termination_condition == TerminationCondition.optimal):
    # Do something when the solution in optimal and feasible
elif (results.solver.termination_condition == TerminationCondition.infeasible):
    # Do something when model in infeasible
else:
    # Something else is wrong
    print “Solver Status: ”,  result.solver.status
'''
