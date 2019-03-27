# -*- coding: utf-8 -*-
from pyomo.environ import *
from pyomo.opt import SolverFactory
import ruizC as RC
import math


############################################################
#SubProblem

#Create an abstract model in Pyomo for the Sub Problem
#Input is expected to be .dat file
#Will take x_star (binary) as a parameter
############################################################

mod = AbstractModel()
opt = SolverFactory(RC.SOLVER)

###############################################################
#Parameters and Sets
###############################################################

#Sets
mod.N = 	Set()						#Nodes
mod.L = 	Set(within=mod.N*mod.N)		#Lines

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

#B Matrix
#Split into parts for each variable
mod.gen_mu = 	Param(mod.N*mod.constraints)
mod.alpha_mu = 	Param(mod.N*mod.constraints)
mod.tran_mu = 	Param(RangeSet(1,mod.LLen)*mod.constraints)

#Parameters that come from Master
mod.x_star = 	Param(mod.L, domain=NonNegativeIntegers, default=0,
	mutable = True) 					#Built Lines


###############################################################
#Variables
###############################################################

#Variables
mod.tran = Var(mod.L, within=Reals) 			#Ammount Transmitted
mod.gen  =	Var(mod.N, domain=NonNegativeReals)	#Generation
mod.genpos = Var(mod.N, domain=NonNegativeReals)#Max Possible Generation
mod.alpha = Var(mod.N, domain=NonNegativeReals)	#Unfilled Demand
mod.dem  =	Var(mod.N, domain=NonNegativeReals)	#Demand
mod.lamda = Var(mod.constraints)				#Lamda Free
mod.mu	=	Var(mod.constraints, domain=NonNegativeReals)	#Mu >= 0
mod.z	=	Var(mod.constraints, domain=Binary)	#For Slackness

#############################################
#Functions
#############################################

#Objective Function
#	max [gen_cost * generation + shed_cost * unfilled_demand]
def obj_expression(mod):
	return (mod.sigma * (sum(mod.gencost[i] * mod.gen[i] for i in mod.N) 
		 + sum(mod.shed[i] * mod.alpha[i] for i in mod.N))
		 + sum(mod.c[j] * mod.x_star[j] for j in mod.L))
mod.Obj = Objective(rule=obj_expression, sense = maximize)


#Transmisson Capacity
#	abs(transmission) <= capacity * lines
def cap_rule(mod, i, j):
	return mod.tran[i,j] <=  mod.cap[i,j] * mod.x_star[i,j]
mod.CapConstraint = Constraint(mod.L, rule=cap_rule)
def cap_rule2(mod, i, j):
	return mod.tran[i,j] >=  -mod.cap[i,j] * mod.x_star[i,j]
mod.CapConstraint2 = Constraint(mod.L, rule=cap_rule2)
#	<= z*M
def cap_rule_zm1(mod, i, j):
	zmem = mod.LLen
	n = mod.NLen
	return (- mod.tran[i,j] + mod.cap[i,j] * mod.x_star[i,j] <=
		mod.M * mod.z[zmem - (n-i)*(n-i-1)/2 + j - i - 1])
mod.CapConstraintZM1 = Constraint(mod.L, rule=cap_rule_zm1)
def cap_rule_zm2(mod, i, j):
	zmem = 2*mod.LLen
	n = mod.NLen
	return (mod.tran[i,j] + mod.cap[i,j] * mod.x_star[i,j] <=
		mod.M *  mod.z[zmem - (n-i)*(n-i-1)/2 + j - i - 1])
mod.CapConstraintZM2= Constraint(mod.L, rule=cap_rule_zm2)



#Unmet Demand is less than Demand
#	alpha("unmet demand") <= demand
def alpha_rule(mod,i):
	return mod.alpha[i] <= mod.dem[i]
mod.AlphaConstraint = Constraint(mod.N, rule=alpha_rule)
#	<= ZM
def alpha_ruleZM(mod,i):
	zmem = 4*mod.LLen + mod.NLen
	return mod.dem[i] - mod.alpha[i] <= mod.M * mod.z[zmem + i]
mod.AlphaConstraintZM2 = Constraint(mod.N, rule=alpha_ruleZM)


#Flow (Supply and Demand)
#Alpha is unfilled demand
#	flow + supply - demand >= alpha
#flow is positive if from Node 1 to 2, and negative if Node 2 to 1
#	alpha <= demand	
def flow_rule(mod, i):
	flowcol = sum(mod.tran[j,j2] for (j,j2) in mod.L
		if (j2 == i)) 
	flowrow = sum(mod.tran[j,j2] for (j,j2) in mod.L
		if (j == i))
	return (flowcol + mod.gen[i] - flowrow - mod.dem[i]
		>=  -mod.alpha[i])
mod.FlowConstraint = Constraint(mod.N, rule=flow_rule)
#	<=ZM
def flow_ruleZM(mod, i):
	zmem = 4*mod.LLen + 2*mod.NLen
	flowcol = sum(mod.tran[j,j2] for (j,j2) in mod.L
		if (j2 == i)) 
	flowrow = sum(mod.tran[j,j2] for (j,j2) in mod.L
		if (j == i))
	return (flowcol + mod.gen[i] - flowrow - mod.dem[i] +  mod.alpha[i]
		<=  mod.M * mod.z[zmem + i])
mod.FlowConstraintZM = Constraint(mod.N, rule=flow_ruleZM)

#Supply min and max
#	generation <= possible_gen <= gen_max
def genpos_rule(mod, i):
	return (mod.genpos[i] <= mod.supmax[i])
mod.GenPosConstraint = Constraint(mod.N, rule=genpos_rule)
def max_sup_rule(mod, i):
	return (mod.gen[i] <= mod.genpos[i])
mod.MaxSupConstraint = Constraint(mod.N, rule=max_sup_rule)

#def min_sup_rule(mod, i):
#	return (mod.supmin[i] <= mod.gen[i])
#mod.MinSupConstraint = Constraint(mod.N, rule=min_sup_rule)


#Demand Min and Max
#	demand_min <= demand <= demand_max
def max_dem_rule(mod, i):
	return mod.dem[i]  <= mod.demmax[i]
mod.MaxDemConstraint = Constraint(mod.N, rule=max_dem_rule)
def min_dem_rule(mod, i):
	return mod.dem[i]  >= mod.demmin[i]
mod.MinDemConstraint = Constraint(mod.N, rule=min_dem_rule)


#Adjust the supply uncertainty budget
def unc_sup_rule(mod):
	if sum(mod.supmax[i] - mod.supmin[i] for i in mod.N) <= RC.TOL:
		return Constraint.Feasible
	else:
		return  (sum(mod.supmax[i] - mod.genpos[i] for i in mod.N)
			/ sum(mod.supmax[i] - mod.supmin[i] for i in mod.N)
			== mod.uncS)
mod.UncSupConstraint = Constraint(rule=unc_sup_rule)


#Adjust the demand uncertainty budget
def unc_dem_rule(mod):
	if sum(mod.demmax[i] - mod.demmin[i] for i in mod.N) <= RC.TOL:
		return Constraint.Feasible
	else:
		return  (sum(mod.dem[i] - mod.demmin[i] for i in mod.N)
			/ sum(mod.demmax[i] - mod.demmin[i] for i in mod.N)
			>= mod.uncD)
mod.UncDemConstraint = Constraint(rule=unc_dem_rule)


#Mu is 0 when z is 1	
def mu_bigm_rule(mod, i):
	return mod.mu[i] <= mod.M*(1-mod.z[i]) 
mod.MuBigM = Constraint(mod.constraints, rule=mu_bigm_rule)

#Lambda Mu Yp Rules
#############################################
#Generation
#I = 6
def lamda_mu_gen_rule(mod, i):
	return (sum(mod.gen_mu[i,j] * mod.mu[j] for j in mod.constraints)
		== -mod.gencost[i])
mod.LamdaMuGenConstraint = Constraint(mod.N, rule=lamda_mu_gen_rule)

#Alpha
# I = 6
def lamda_mu_alpha_rule(mod, i):
	return (sum(mod.alpha_mu[i,j] * mod.mu[j] for j in mod.constraints)
		== -mod.shed[i])
mod.LamdaMuAlphaConstraint = Constraint(mod.N, rule=lamda_mu_alpha_rule)

#Transmission
# I = 15
def lamda_mu_tran_rule(mod, i):
	return (sum(mod.tran_mu[i,j] * mod.mu[j] for j in mod.constraints)
		== 0)
mod.LamdaMuTranConstraint = Constraint(mod.LRange, rule=lamda_mu_tran_rule)
#############################################




 

##########
#TO TEST
##########

'''
isub = mod.create_instance(RC.DATA)
results = opt.solve(isub, tee=True)
#isub.pprint()
results.write()

##To Print

for v in isub.component_objects(Var, active=True):
	print ("Variable",v)
	varob = getattr(isub, str(v))
	for index in varob:
		print ("   ",index, varob[index].value) 
'''

