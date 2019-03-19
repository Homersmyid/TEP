# -*- coding: utf-8 -*-
from pyomo.environ import *
from pyomo.opt import SolverFactory
import ruizC as RC

###############################################################
#Master Problem Abstract Model

#Create an abstract model in Pyomo for the Master Problem
###############################################################

mod = AbstractModel()					#name of model
opt = SolverFactory(RC.SOLVER)			#solver to use

###############################################################
#Parameters and Variables
###############################################################
#Parameters
mod.N = 	Set()						#Nodes
mod.L = 	Set(within=mod.N*mod.N)		#Lines
mod.c = 	Param(mod.L)				#cost per line
mod.pi = Param()						#Budget
mod.maxLines = Param()					#Max Lines per route
mod.cap =	Param(mod.L)				#line capacity
mod.sigma = Param()						#Hours in a year
mod.demmax = Param(mod.N)				#Maximum possible Demand
mod.demmin = Param(mod.N)				#Minimum possible Demand
mod.supmax = Param(mod.N)				#Maximum possible Supply
mod.supmin = Param(mod.N)				#Minimum possible Supply
mod.M	=	Param()						#Max M for Fourtuny-Amat
mod.gencost =	Param(mod.N)			#Cost to generate	
mod.shed =  Param(mod.N)				#Load Shedding Cost Per Node
mod.uncD =  Param()						#Uncertainty in Demand	
mod.uncS =  Param()						#Uncertainty in Supply
mod.conLen = Param()					#Constraints in Primal 
mod.constraints = RangeSet(1,mod.conLen) 	#(1, '# of constraints')
mod.varLen = Param()					#Variables in Primal
mod.NLen = Param()						#Length of N
mod.LLen = Param()						#Length of L
mod.LRange = RangeSet(1,mod.LLen)		#(1, '# of lines')

#B Matrix
#Split into parts for each variable
mod.gen_mu = Param(mod.N*mod.constraints)
mod.alpha_mu = Param(mod.N*mod.constraints)
mod.tran_mu = Param(RangeSet(1,mod.LLen)*mod.constraints)

#Parameters that come from subproblem
mod.dem = 	 Param(mod.N, default=0, mutable = True)	#Demand
mod.genpos = Param(mod.N, default=0, mutable = True)	#Maximum Supply
mod.x_star = Param(mod.L, domain=NonNegativeIntegers, default=0,
	mutable = True) 									#Built Lines

#Variables
mod.x 	 = 	Var(mod.L, domain=NonNegativeIntegers)	#Lines Built
mod.eta = 	Var(domain=NonNegativeReals)	#Eta >= b^t*(yp) for all p
mod.tran = Var(mod.L, within=Reals) 				#Ammount Transmitted
mod.gen  =	Var(mod.N, domain=NonNegativeReals)		#Generation Supply
mod.alpha = Var(mod.N, domain=NonNegativeReals)		#Unfilled Demand

###############################################################
#Functions
###############################################################

#Objective Function
# 	min [c^t*x + sigma * eta]
def obj_expression(mod):
	return (sum(mod.c[i,j] * mod.x[i,j] for i,j in mod.L)
		+ mod.sigma * mod.eta)
mod.Obj = Objective(rule=obj_expression)


#Budget Constraint (on power line building)
#	c^t*x <= pi
def budget_constraint_rule(mod):
	return summation(mod.c,mod.x) <= mod.pi
mod.BudgetConstraint = Constraint(rule=budget_constraint_rule)


#X_star
#x has to have more than or eqaul to lines that x_star has
def x_star_rule(mod, i, j):
	return  mod.x[i, j] >= mod.x_star[i, j]
mod.XStarConstraint = Constraint(mod.L, rule=x_star_rule)


#Max Lines per connection
def max_lines_rule(mod, i, j):
	return  mod.x[i, j] <= mod.maxLines
mod.MaxLinesRuleConstraint = Constraint(mod.L, rule=max_lines_rule)


#Transmisson Capacity
#	abs(transmission) <= capacity
def cap_rule(mod, i, j):
	return mod.tran[i,j] <=  mod.cap[i,j] * mod.x[i,j]
mod.CapConstraint = Constraint(mod.L, rule=cap_rule)
def cap_rule2(mod, i, j):
	return -mod.tran[i,j] <=  mod.cap[i,j]* mod.x[i,j]
mod.CapConstraint2 = Constraint(mod.L, rule=cap_rule2)

#Supply min and max
#	Yp <= possible_generation
def gen_rule(mod, i):
	return (mod.gen[i] <= mod.genpos[i])
mod.GenConstraint = Constraint(mod.N, rule=gen_rule)


#Unmet Demand is less than Demand
#	alpha("unmet demand") <= demand
def alpha_rule(mod,i):
	return mod.alpha[i] <= mod.dem[i]
mod.AlphaConstraint = Constraint(mod.N, rule=alpha_rule)


#Flow (Supply and Demand)
#	flow + supply - demand >= alpha
#Alpha is unmet demand
#Flow is positive if from Node 1 to 2, and negative if Node 2 to 1
def flow_rule(mod, i):
	flowcol = sum(mod.tran[j,j2] for (j,j2) in mod.L
		if (j2 == i)) 
	flowrow = sum(mod.tran[j,j2] for (j,j2) in mod.L
		if (j == i))
	return (flowcol + mod.gen[i] - flowrow - mod.dem[i]
		>=  -mod.alpha[i])
mod.FlowConstraint = Constraint(mod.N, rule=flow_rule)


# Hourly Costs 
#	Eta >= b^t*(y) for all yp
#Costs are :
#	1) generation costs
#	2) penalty for unmet demand
def eta_rule(mod):
	return (sum(mod.gencost[i] * mod.gen[i] for i in mod.N) 
		 +sum(mod.shed[i] * mod.alpha[i] for i in mod.N) <= mod.eta)
mod.EtaConstraint = Constraint(rule=eta_rule)




##########
#TO TEST
##########

'''
imast = mod.create_instance(RC.DATA)	
results = opt.solve(imast, tee=True)
#imast.pprint()
results.write()


#To Print
for v in imast.component_objects(Var, active=True):
	print ("Variable",v)
	varob = getattr(imast, str(v))
	for index in varob:
		print ("   ",index, varob[index].value) 

'''
