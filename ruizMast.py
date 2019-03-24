# -*- coding: utf-8 -*-
from pyomo.environ import *
from pyomo.opt import SolverFactory
import ruizC as RC
import math

###############################################################
#Master Problem Abstract Model
#Create an abstract model in Pyomo for the Master Problem
###############################################################

mod = AbstractModel()					#name of model
opt = SolverFactory(RC.SOLVER)			#solver to use
mod.con = ConstraintList()

###############################################################
#Parameters and Sets
###############################################################

#Sets
mod.N = 	Set()						#Nodes
mod.L = 	Set(within=mod.N*mod.N)		#Lines
mod.P = 	Set(initialize=[])			#Subproblem Solves
	#Updated every new run of the master

#Parameters
mod.c = 	Param(mod.L)				#Cost per line
mod.pi = Param()						#Budget
mod.maxLines = Param()					#Max Lines per route
mod.cap =	Param(mod.L)				#Line capacity
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
mod.b =		Param(mod.L)				#Phyiscs on each line
mod.Mtheta = Param()					#M to use for theta constraint
mod.ref	=	Param(mod.N)				#Reference Theta

#Parameters that come from subproblem
mod.dem = 	 Param(mod.N, default=0, mutable = True)	#Demand
mod.genpos = Param(mod.N, default=0, mutable = True)	#Maximum Supply
mod.x_star = Param(mod.L, domain=NonNegativeIntegers, default=0,
	mutable = True) 									#Built Lines

###############################################################
#Variables
###############################################################

#Variables
mod.x 	 = 	Var(mod.L, domain=NonNegativeIntegers)	#Lines Built
mod.eta = 	Var(domain=NonNegativeReals)	#Eta >= b^t*(yp) for all p
mod.tran = Var(mod. P, mod.L, within=Reals) 		#Ammount Transmitted
mod.gen  =	Var(mod.N, domain=NonNegativeReals)		#Generation Supply
mod.alpha = Var(mod.N, domain=NonNegativeReals)		#Unfilled Demand
mod.theta = Var(mod.N,bounds=(-math.pi, math.pi))	#Angle in [-pi,pi]
mod.route_on =	Var(mod.L, domain=Binary)			#If route used 
	#Needed for when number of lines per route is >= 1


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


# Route Rule
# To see what routes are activated, route_on is binary
# MaxLines * Route_on >= x
def route_rule(mod, i,j):
	return (mod.maxLines * mod.route_on[i,j] >= mod.x[i,j])
mod.RouteConstraint = Constraint(mod.L, rule=route_rule)









'''
#Transmisson Capacity
#	abs(transmission) <= capacity
def cap_rule1(mod, i, j):
	return mod.tran[i,j] <=  mod.cap[i,j] * mod.x[i,j]
mod.CapConstraint1 = Constraint(mod.L, rule=cap_rule1)
def cap_rule2(mod, i, j):
	return -mod.tran[i,j] <=  mod.cap[i,j]* mod.x[i,j]
mod.CapConstraint2 = Constraint(mod.L, rule=cap_rule2)


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


# Eta Rule (Hourly Costs)
#	Eta >= b^t*(y) for all yp
#Costs are :
#	1) generation costs
#	2) penalty for unmet demand
def eta_rule(mod):
	return (sum(mod.gencost[i] * mod.gen[i] for i in mod.N) 
		 +sum(mod.shed[i] * mod.alpha[i] for i in mod.N) <= mod.eta)
mod.EtaConstraint = Constraint(rule=eta_rule)


# Theta Rules
def theta_rule1(mod,i,j):
	return (mod.tran[i,j] - mod.b[i,j] * (mod.theta[i] - mod.theta[j])
		<= (1 - mod.route_on[i,j]) * mod.Mtheta)
mod.ThetaConstraint1 = Constraint(mod.L, rule=theta_rule1)
def theta_rule2(mod,i,j):
	return (-(mod.tran[i,j] - mod.b[i,j] * (mod.theta[i] - mod.theta[j]))
		<= (1 - mod.route_on[i,j]) * mod.Mtheta)
mod.ThetaConstraint2 = Constraint(mod.L, rule=theta_rule2)

######################
#NEEDED? will this just slow donw program

# Set a reference theta to zero 
# |theta| <= ref
def ref_rule1(mod,i):
	return (mod.theta[i] <= mod.ref[i]) 
mod.RefConstraint1 = Constraint(mod.N, rule=ref_rule1)
def ref_rule2(mod,i):
	return (-mod.theta[i] <= mod.ref[i]) 
mod.RefConstraint2 = Constraint(mod.N, rule=ref_rule2)
#######################
'''

##########
#TO TEST
##########

'''
imast = mod.create_instance(RC.DATA)

#Set x_star in step zero
for x in RC.START_X_STAR:
	imast.x_star[x[0], x[1]] = x[2]

#set demand in master
list1 = [200, 0,0, 150, 100, 200]
for d in range(1,7):
	imast.dem[d] = list1[d-1]

#set possible generation in master
list2 = [300, 250, 400, 0, 300, 150]
for d in range(1,7):
	imast.genpos[d] = list2[d-1]
	
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

###############################################################
#Functions
###############################################################

def mast_func(imast, subdem, subgenpos):

	#Supply min and max
	#	Gen <= (Possible_Generation)
	for i in subgenpos:
		imast.con.add(imast.gen[i] <= value(subgenpos[i]))
		
	#Unmet Demand is less than Demand
	#	Alpha("unmet demand") <= demand
	for i in subdem:
		imast.con.add(imast.alpha[i] <= value(subdem[i]))
	

