# -*- coding: utf-8 -*-
from pyomo.environ import *
from pyomo.opt import SolverFactory
import ruizC as RC


###############################################################
#Master Problem Abstract Model

#Create an abstract model in Pyomo for the Master Problem
#Input is expected to be .dat file
#Variables x are assumed to be bianry
###############################################################

mast = AbstractModel()
opt = SolverFactory(RC.SOLVER)

#Parameters
mast.N = 	Set()							#Nodes
mast.L = 	Set(within=mast.N*mast.N)		#Lines
mast.c = 	Param(mast.L)					#cost per line
mast.pi = 	Param()							#budget for power lines
mast.sigma = Param()						#Hours in a year (IF NEEDED)
mast.cap =	Param(mast.L)					#line capacity
mast.gencost =	Param(mast.N)				#Cost to generate	
mast.shed =  Param(mast.N)					#Load Shedding Cost Per Node
mast.maxLines = Param()						#Max Lines Per Connection

#Parameters that come from subproblem
mast.dem = Param(mast.N, default=0, mutable = True)		#Demand
mast.genpos = Param(mast.N, default=0, mutable = True)	#Maximum Supply
mast.x_star = Param(mast.L, domain=NonNegativeIntegers, default=0,
	mutable = True) 									#x_star

#Variables
mast.x 	 = 	Var(mast.L, domain=NonNegativeIntegers)	#Lines Built
mast.eta = 	Var(domain=NonNegativeReals)	# Eta >= b^t*(yp) for all p
mast.tran = Var(mast.L, within=Reals) 		#Ammount Transmitted
mast.gen  =	Var(mast.N, domain=NonNegativeReals)	#Generation Supply
mast.alpha = Var(mast.N, domain=NonNegativeReals)	#Unfilled Demand

###############################################################
#Functions
###############################################################

#Objective Function
# 	min [c^t*x + sigma * eta]
def obj_expression(mast):
	return (sum(mast.c[i,j] * mast.x[i,j] for i,j in mast.L)
		+ mast.sigma * mast.eta)
mast.OBJ = Objective(rule=obj_expression)


#Budget Constraint (on power line building)
#	c^t*x <= pi
def budget_constraint_rule(mast):
	return summation(mast.c,mast.x) <= mast.pi
mast.BudgetConstraint = Constraint(rule=budget_constraint_rule)


#X_star
#x has to have more than or eqaul to lines that x_star has
def x_star_rule(mast, i, j):
	return  mast.x[i, j] >= mast.x_star[i, j]
mast.XStarConstraint = Constraint(mast.L, rule=x_star_rule)

#Max Lines per connection
def max_lines_rule(mast, i, j):
	return  mast.x[i, j] <= mast.maxLines
mast.MaxLinesRuleConstraint = Constraint(mast.L, rule=max_lines_rule)

#Transmisson Capacity
#	abs(transmission) <= capacity
def cap_rule(mast, i, j):
	return mast.tran[i,j] <=  mast.cap[i,j] * mast.x[i,j]
mast.CapConstraint = Constraint(mast.L, rule=cap_rule)
def cap_rule2(mast, i, j):
	return -mast.tran[i,j] <=  mast.cap[i,j]* mast.x[i,j]
mast.CapConstraint2 = Constraint(mast.L, rule=cap_rule2)

#Supply min and max
#	Yp <= possible_generation
def gen_rule(mast, i):
	return (mast.gen[i] <= mast.genpos[i])
mast.GenConstraint = Constraint(mast.N, rule=gen_rule)


#Unmet Demand is less than Demand
#	alpha("unmet demand") <= demand
def alpha_rule(mast,i):
	return mast.alpha[i] <= mast.dem[i]
mast.AlphaConstraint = Constraint(mast.N, rule=alpha_rule)


#Flow (Supply and Demand)
#	flow + supply - demand >= alpha
#Alpha is unmet demand
#Flow is positive if from Node 1 to 2, and negative if Node 2 to 1
def flow_rule(mast, i):
	flowcol = sum(mast.tran[j,j2] for (j,j2) in mast.L
		if (j2 == i)) 
	flowrow = sum(mast.tran[j,j2] for (j,j2) in mast.L
		if (j == i))
	return (flowcol + mast.gen[i] - flowrow - mast.dem[i]
		>=  -mast.alpha[i])
mast.FlowConstraint = Constraint(mast.N, rule=flow_rule)


# Hourly Costs are kept are seperate "eta" for convience
#	Eta >= b^t*(y)
#Costs are :
#	1) generation costs
#	2) penalty for unmet demand
def eta_rule(mast):
	return (sum(mast.gencost[i] * mast.gen[i] for i in mast.N) 
		 +sum(mast.shed[i] * mast.alpha[i] for i in mast.N) <= mast.eta)
mast.EtaConstraint = Constraint(rule=eta_rule)




##########
#TO TEST
##########
'''
imast = mast.create_instance(RC.MAST)	
results = opt.solve(imast, tee=True)
imast.pprint()
results.write()

#To Print
for v in imast.component_objects(Var, active=True):
	print ("Variable",v)
	varob = getattr(imast, str(v))
	for index in varob:
		print ("   ",index, varob[index].value) 

'''
