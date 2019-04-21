# -*- coding: utf-8 -*-

from __future__ import division
from pyomo.environ import *
from pyomo.opt import SolverFactory
import ruizC as RC
import ruizSub as s
import ruizMast as m

STOP = 4;							#How many iterations to quit after
startlines = True					#If possible lines at start
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
	imast = m.mod.create_instance(RC.DATA)
		
	#Set x_star in step zero
	for x in RC.START_X_STAR:
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
	if abs(UB - LB) / UB <= RC.EPSILON:
		break
	
	########################
	#STEP K Master Problem
	######################## 
	for x in [(1,300),(2,250),(3,400),(4,0),(5,170),(6,60)]:
		isub.gen[x[0]] = x[1]
	m.mast_func(imast, isub.dem,isub.gen, RC.START_X_STAR, k)

	#solve master problem
	mresults = s.opt.solve(imast)
	LB = value(imast.Obj)
	
	
	
	print('\n\nk:', k)
	print("*MASTER*\n\n")
	mresults.write()	
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
	
	#store new upper bound if it is < previous upper bound
	if value (isub.Obj) <= UB:
		UB = value(isub.Obj)
	
	#isub.pprint()
	#input()

	print('\n\nk:', k)
	print("*SUB***\n\n")
	sresults.write()
	for v in isub.component_objects(Var, active=True):
		print ("Variable",v)
		varob = getattr(isub, str(v))
		for index in varob:
			print ("   ",index, varob[index].value)
	input()
	
	print("XXX")
	print(UB)
	print(LB)
	print("XXX")

