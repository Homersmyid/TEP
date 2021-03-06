######################################################################################################################################
#Implements the master problem of the column-and-constraint method.
#Has two components, the abstract model in Pyomo of the master problem
#and a function for adding variables and constraints to the concrete
#model on every new iteration.
#See file MastTheory for a full description of the theory behind setup.
#######################################################################


#######################################################################
#Master Problem Abstract Model

#This method takes in the realizations of uncertain parameters, both
#supply and demand, from previous subproblems as constant parameters.
#Then, it outputs the choice of lines that will minimize both line
#construction costs and hourly costs, as well a new lower bound,
#which will always be less than or equal to previous lower bound.

#Input–	(From user) None
#		(From main) X_star. A list of already existing power lines.
#       (From data file) Data for concrete model. See docs/Input.pdf.

#Output- (To main) Results of a Pyomo.opt.solve. 
#			       Value of maximization solve as a lower bound.
#######################################################################
#Master Problem Abstract Model
#Create an abstract model in Pyomo for the Master Problem
###############################################################

# -*- coding: utf-8 -*-
from pyomo.environ import *
from pyomo.opt import SolverFactory
import ruizC as RC
import math

mod = AbstractModel()					#name of model
opt = SolverFactory(RC.SOLVER)			#solver to use
opt.options['mipgap'] = RC.MIPGAP

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
mod.cap =		Param(mod.L)			#Line capacity
mod.b =			Param(mod.L)			#Phyiscs on each line
mod.pi = 		Param()					#Budget
mod.maxLines =	Param()					#Max Lines per route
mod.sigma = 	Param()					#Hours in a year
mod.demmax = 	Param(mod.N)			#Maximum possible Demand
mod.demmin = 	Param(mod.N)			#Minimum possible Demand
mod.supmax = 	Param(mod.N)			#Maximum possible Supply
mod.supmin = 	Param(mod.N)			#Minimum possible Supply
mod.gencost =	Param(mod.N)			#Cost to generate	
mod.shed =  	Param(mod.N)			#Load Shedding Cost Per Node
mod.uncD =  	Param()					#Uncertainty in Demand	
mod.uncS =  	Param()					#Uncertainty in Supply
mod.ref	=		Param()					#Reference Theta

#Big Ms
mod.M	=		Param()					#Max M for Duals
mod.Mgen =		Param()					#Highest Gen Possible
mod.Mdem =		Param()					#Highest Demand Possible
mod.Mcap = 		Param()					#Highest Line Capacity
mod.Mtheta = 	Param()					#7 since 7 > 2pi
mod.Mtran = 	Param()					#Most that can be transmitted
	#Note if line more than 1 than this can be more than any capacity

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

#############################
	
#Variables that grow each subproblem solve (P = subproblem solves)

#Ammount Transmitted
mod.tran = Var(mod.P, mod.L, within=Reals, initialize=0)

#Generation Supply	
mod.gen  =	Var(mod.P, mod.N, domain=NonNegativeReals, initialize=0)
	
#Unfilled Demand	
mod.unmet = Var(mod.P, mod.N, domain=NonNegativeReals, initialize=0)

#Angle in [-pi,pi]	
mod.theta = Var(mod.P, mod.N, bounds=(-math.pi, math.pi), initialize=0)	

###############################################################
#Constraints
###############################################################

#Objective Function
# 	min [c^t*x + eta]
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
# 	X <= MaxLines * Route_on <= MaxLines * X 
def route_rule1(mod, i,j):
	return (mod.x[i,j] <= mod.maxLines * mod.route_on[i,j])
mod.RouteConstraint1 = Constraint(mod.L, rule=route_rule1)
def route_rule2(mod, i,j):
	return (mod.route_on[i,j] <=  mod.x[i,j])
mod.RouteConstraint2 = Constraint(mod.L, rule=route_rule2)

###############################################################
#Expanding Constraints
# A copy of these will be added for each new run of the master
###############################################################
mod.GenConstraint 		= ConstraintList()
mod.UnmetDemConstraint 	= ConstraintList()
mod.CapConstraintPos 	= ConstraintList()
mod.CapConstraintNeg 	= ConstraintList()
mod.FlowConstraint 		= ConstraintList()
mod.ThetaConstraint 	= ConstraintList()
mod.EtaConstraint 		= ConstraintList()
mod.RefConstraint 		= ConstraintList()





###############################################################
#Master Function

# mast_func(imast, subdem, subgenpos, in_x_star, k)

#Input- Imast. The concrete version of the master problem.
#		Subdem. Demand at each sink. From last subproblem solve.
#		Subgenpos. Generation capacity at each generator. From last subproblem solve.
#		In_x_star. A list of preexisting power lines.
#		K. What iteration the main is currently on.

#Output- Makes changes in the "imast" function
#
# This function adds a new series of constraints based on the
# generation and supply level found from the previous subproblem
# solve. 
###############################################################

def mast_func(imast, subdem, subgenpos, in_x_star, k):
	
	#Add to set P that there is a new subproblem solved
	imast.P.add(k)
	
	#Set demand in master
	for i in subdem:
		imast.dem[i] = value(subdem[i])

	#Set possible generation in master
	for i in subgenpos:
		imast.genpos[i] = value(subgenpos[i])

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
		flowcol = sum(imast.tran[k,j,j2] for (j,j2) in imast.L
			if (j2 == i))
		flowrow = sum(imast.tran[k,j,j2] for (j,j2) in imast.L
			if (j == i))
		imast.FlowConstraint.add( flowcol - flowrow
			+ imast.gen[k,i] - imast.dem[i] >= -imast.unmet[k,i])

	# Theta Rules
	# Theta is the angle at each node
	# For each route from i to j:
	#	(x)*(b)*(theta_i - theta_j) = flow  (if route "ij" active)		
	for i,j in imast.L:
		imast.ThetaConstraint.add( imast.x_star[i,j] * imast.b[i,j]
		* (imast.theta[k,i] - imast.theta[k,j]) == imast.tran[k,i,j])
	
	# Eta Rule (Hourly Costs)
	#	Eta >= sigma * [(Generation_Costs) + (Load_Shedding)]
	# Sigma is hours in a year
	imast.EtaConstraint.add( imast.sigma * 
		(sum(imast.gencost[i] * imast.gen[k,i] for i in imast.N) 
		+ sum(imast.shed[i] * imast.unmet[k,i] for i in imast.N))
		<= imast.eta)

	# Reference Theta
	#	Theta refernce = 0 for each k 
	imast.RefConstraint.add(imast.theta[k,value(imast.ref)] == 0)
	

