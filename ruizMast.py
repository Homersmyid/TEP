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

###############################################################
#Parameters and Sets
###############################################################

#Sets
mod.N = 	Set()						#Nodes
mod.L = 	Set(within=mod.N*mod.N)		#Lines
mod.P = 	Set(initialize=[])			#Subproblem Solves
	#Updated every new run of the master

#Parameters
mod.c = 		Param(mod.L)			#Cost per line
mod.pi = 		Param()					#Budget
mod.maxLines =	Param()					#Max Lines per route
mod.cap =		Param(mod.L)			#Line capacity
mod.sigma = 	Param()					#Hours in a year
mod.demmax = 	Param(mod.N)			#Maximum possible Demand
mod.demmin = 	Param(mod.N)			#Minimum possible Demand
mod.supmax = 	Param(mod.N)			#Maximum possible Supply
mod.supmin = 	Param(mod.N)			#Minimum possible Supply
mod.M	=		Param()					#Max M for Fourtuny-Amat
mod.gencost =	Param(mod.N)			#Cost to generate	
mod.shed =  	Param(mod.N)			#Load Shedding Cost Per Node
mod.uncD =  	Param()					#Uncertainty in Demand	
mod.uncS =  	Param()					#Uncertainty in Supply
mod.conLen = 	Param()					#Constraints in Primal 
mod.constraints = RangeSet(1,mod.conLen) 	#(1, '# of constraints')
mod.varLen = 	Param()					#Variables in Primal
mod.NLen = 		Param()					#Length of N
mod.LLen = 		Param()					#Length of L
mod.LRange = RangeSet(1,mod.LLen)		#(1, '# of lines')
mod.b =			Param(mod.L)			#Phyiscs on each line
mod.Mtheta = 	Param()					#M to use for theta constraint
mod.ref	=		Param(mod.N)			#Reference Theta

#Parameters that come from subproblem or setup
mod.dem = 	 Param(mod.N, default=0, mutable = True)	#Demand
mod.genpos = Param(mod.N, default=0, mutable = True)	#Maximum Supply
mod.x_star = Param(mod.L, domain=NonNegativeIntegers, default=0,
	mutable = True) 									#Built Lines

###############################################################
#Variables
###############################################################

#Variables

#Lines Built on each route (Integer)
mod.x 	 = 	Var(mod.L, domain=NonNegativeIntegers, initialize=0)

#Worst Case Hourly Costs (Generation + Load Shed)
mod.eta = 	Var(domain=NonNegativeReals, initialize=0)

#If route used (Bianary). Needed for when the lines per route is > 1
mod.route_on =	Var(mod.L, domain=Binary, initialize=0)			

	
#Variables that grow each subproblem solve (P = subproblem solves)

#Ammount Transmitted
mod.tran = Var(mod. P, mod.L, within=Reals, initialize=0)

#Generation Supply	
mod.gen  =	Var(mod. P, mod.N, domain=NonNegativeReals, initialize=0)
	
#Unfilled Demand	
mod.unmet = Var(mod. P, mod.N, domain=NonNegativeReals, initialize=0)

#Angle in [-pi,pi]	
mod.theta = Var(mod. P, mod.N,bounds=(-math.pi, math.pi), initialize=0)	

###############################################################
#Constraints
###############################################################

#Objective Function
# 	min [c^t*x + sigma * eta]
def obj_expression(mod):
	return (sum(mod.c[i,j] * mod.x[i,j] for i,j in mod.L)
		+ mod.eta)
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


###############################################################
#Expanding Constraints

#A copy of these will be added for each new run of the master
###############################################################

mod.GenConstraint = ConstraintList()
mod.UnmetDemConstraint = ConstraintList()
mod.CapConstraintPos = ConstraintList()
mod.CapConstraintNeg = ConstraintList()
mod.FlowConstraint = ConstraintList()
mod.ThetaConstraintPos = ConstraintList()
mod.ThetaConstraintNeg = ConstraintList()
mod.EtaConstraint = ConstraintList()






###############################################################
#Function
###############################################################

def mast_func(imast, subdem, subgenpos, in_x_star, k):

	#Add to set P that there is a new subproblem solved
	imast.P.add(k)
	
	print("...")
	#Set demand in master
	for i in subdem:
		imast.dem[i] = value(subdem[i])
		print(imast.dem[i].value)

	#Set possible generation in master
	for i in subgenpos:
		imast.genpos[i] = value(subgenpos[i])
		
		
		print('a')
		print(subgenpos[i].value)
		print('b')	
		print(imast.genpos[i].value)
	print("***")


	#Set x_star
	for x in in_x_star:
		imast.x_star[x[0], x[1]] = x[2]

	#Supply min and max
	#	Gen <= (Possible_Generation)
	for i in subgenpos:
		imast.GenConstraint.add( imast.gen[k,i] <= value(subgenpos[i]))
		
	#Unmet Demand is less than Demand
	#	("Unmet Demand") <= demand
	for i in subdem:
		imast.UnmetDemConstraint.add( imast.unmet[k,i]
			<= value(subdem[i]))

	#Transmisson Capacity
	#	abs(transmission) <= capacity
	for i in imast.cap:
		imast.CapConstraintPos.add( imast.tran[k,i]
			<=  imast.cap[i] * imast.x[i])			
		imast.CapConstraintNeg.add( -imast.tran[k,i]
		<=  imast.cap[i] * imast.x[i])


	#Flow (Supply and Demand)
	#	Flow + Generation - Demand >= (-1)*("Unmet Demand")
	#Flow is positive if from Node 1 to 2, and negative if Node 2 to 1	
	for i in imast.N:
		flowcol = 0;
		flowcol = sum(imast.tran[k,j,j2] for (j,j2) in imast.L
			if (j2 == i))
		flowrow = sum(imast.tran[k,j,j2] for (j,j2) in imast.L
			if (j == i))
		imast.FlowConstraint.add( flowcol - flowrow
			+ imast.gen[k,i] - imast.dem[i] >= -imast.unmet[k,i])
 
	# Theta Rules
	# Theta is the angle at each node
	# For each line from i to j:
	#	|Transmision - (b)*(theta_i - theta_j)| <= M(1-route_on)
	for i,j in imast.L:
		imast.ThetaConstraintPos.add( imast.tran[k,i,j]
			- imast.b[i,j] * (imast.theta[k, i] - imast.theta[k, j])
			<= (1 - imast.route_on[i,j]) * imast.Mtheta)
		imast.ThetaConstraintNeg.add ( -imast.tran[k,i,j]
			+ imast.b[i,j] * (imast.theta[k, i] - imast.theta[k, j])
			<= (1 - imast.route_on[i,j]) * imast.Mtheta)
 

	# Eta Rule (Hourly Costs)
	#	Eta >= sigma * [(Generation_Costs) + (Load_Shedding)]
	# Sigma is hours in a year
	imast.EtaConstraint.add( imast.sigma * 
		(sum(imast.gencost[i] * imast.gen[k,i] for i in imast.N) 
		+ sum(imast.shed[i] * imast.unmet[k,i] for i in imast.N))
		<= imast.eta)

		
	'''
	######################
	#NEEDED? will this just slow down program

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
