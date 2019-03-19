# -*- coding: utf-8 -*-
from pyomo.environ import *
from pyomo.opt import SolverFactory
import ruizC as RC

###############################################################
#Step Zero

#Simply Find the intital cost of the lines alread built
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
mod.x_star = Param(mod.L, domain=NonNegativeIntegers, default=0,
	mutable = True) 					#Built Lines
mod.conLen = Param()					#Constraints in Primal 
mod.constraints = RangeSet(1,mod.conLen) #(1, '# of constraints')
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
mod.dem = Param(mod.N, default=0, mutable = True)		#Demand
mod.genpos = Param(mod.N, default=0, mutable = True)	#Maximum Supply
mod.x_star = Param(mod.L, domain=NonNegativeIntegers, default=0,
	mutable = True) 

#Variables
mod.x 	 = 	Var(mod.L, domain=NonNegativeIntegers)	#Lines Built


###############################################################
#Functions
###############################################################

#Objective Function
# 	min [c^t*x]
def obj_expression(mod):
	return (sum(mod.c[i,j] * mod.x[i,j] for i,j in mod.L))
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

