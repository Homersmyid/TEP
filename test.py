# -*- coding: utf-8 -*-

from __future__ import division
from pyomo.environ import *
from pyomo.opt import SolverFactory
import numpy as np
import ruizC as RC
import ruizSub as s
import ruizMast as m
import ruizStep0 as szero

STOP = 5;							#How many iterations to quit after
startlines = True					#If possible lines at start
k = 0								#Iterations
LB = float("-inf")					#Upper Bound
UB = float("inf")					#Lower Bound
genlist = []				#List of generations from each subproblem
demlist = []				#List of demands from each subproblem

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
	LB = value(imast.Obj)
	
	'''
	print("\n\n***MASTER ZERO***\n\n")
	zresults.write()
	for v in imast.component_objects(Var, active=True):
		print ("Variable",v)
		varob = getattr(imast, str(v))
		for index in varob:
			print ("   ",index, varob[index].value)
	input()
	'''
	
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

'''
print("\n\n***SUB ZERO***\n\n")
sresults.write()
for v in isub.component_objects(Var, active=True):
	print ("Variable",v)
	varob = getattr(isub, str(v))
	for index in varob:
		print ("   ",index, varob[index].value)
#isub.pprint()
input()
'''

############################
#Main Loop
############################
for k in range(1,STOP+1):
	
	#If UB and LB close enough, quit loop
	#if abs(UB - LB) / UB <= RC.EPSILON:
	#	break
	
	########################
	#STEP K Master Problem
	######################## 
	#create master problem
	if k == 1:			
		imast = m.mod.create_instance(RC.DATA)
	
	#Add to list of demands for each subproblem
	temp = []
	for i in isub.dem:
		print(value(isub.dem[i]))
		temp.append(value(isub.dem[i]))
	demlist.append(temp)
	
	temp = []
	for i in isub.genpos:
		print(value(isub.genpos[i]))
		temp.append(value(isub.dem[i]))
	genlist.append(temp)
		
	print("Blake")
	print(demlist)
	print(genlist)

######################################################
	#set demand in master
	for d in isub.dem:
		imast.dem[d] = value(isub.dem[d])

	#set possible generation in master
	for g in isub.genpos:
		imast.genpos[g] = value(isub.genpos[g])
		
	#Set x_star in master
	for x in RC.START_X_STAR:
		imast.x_star[x[0], x[1]] = x[2]
#####################################################


	m.mast_func(imast)

	#solve master problem
	mresults = s.opt.solve(imast)
	print('\n\nk:', k)
	print("*MASTER*\n\n")
#	mresults.write()
	LB = value(imast.Obj)
 
'''	
	for v in imast.component_objects(Var, active=True):
		print ("Variable",v)
		varob = getattr(imast, str(v))
		for index in varob:
			print ("   ",index, varob[index].value)
	#imast.pprint();
	input()
'''
'''
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
'''

