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
mod.unmet = Var(mod.N, domain=NonNegativeReals)	#Unfilled Demand
mod.dem  =	Var(mod.N, domain=NonNegativeReals)	#Demand
mod.lamda = Var(mod.constraints)				#Lamda Free
mod.mu	=	Var(mod.constraints, domain=NonNegativeReals)	#Mu >= 0

#Binary Variables for Complementary Condition Linearization
mod.z_unmet	=	Var(mod.N, domain=Binary)	#For Unmet Demand

#Angle in [-pi,pi]
mod.theta = Var(mod.N,bounds=(-math.pi, math.pi), initialize=0)	

#############################################
#Functions
#############################################

#Objective Function
#	max [gen_cost * generation + shed_cost * unfilled_demand]
def obj_expression(mod):
	return (mod.sigma * (sum(mod.gencost[i] * mod.gen[i] for i in mod.N) 
		 + sum(mod.shed[i] * mod.unmet[i] for i in mod.N))
		 + sum(mod.c[j] * mod.x_star[j] for j in mod.L))
mod.Obj = Objective(rule=obj_expression, sense = maximize)


#GenPos min and max
# Gen_pos is generation possible at a node
#	(Generation) <= possible_gen <= gen_max
def max_genpos_rule(mod, i):
	return (mod.genpos[i] <= mod.supmax[i])
mod.MaxGenPosConstraint = Constraint(mod.N, rule=max_genpos_rule)
def min_genpos_rule(mod, i):
	return (mod.genpos[i] >= mod.supmin[i])
mod.MinGenPosConstraint = Constraint(mod.N, rule=min_genpos_rule)

	
#Generation Max
# Genation is the decision variable being minimized by the agent
# Generation <= Generation Possible
def max_gen_rule(mod, i):
	return (mod.gen[i] <= mod.genpos[i])
mod.MaxSupConstraint = Constraint(mod.N, rule=max_gen_rule)

#Generation Dual	[Phi^(E.min), Phi^(E.max)]
#
#
#


#Demand Min and Max
#	demand_min <= (Demand) <= demand_max
def max_dem_rule(mod, i):
	return mod.dem[i]  <= mod.demmax[i]
mod.MaxDemConstraint = Constraint(mod.N, rule=max_dem_rule)
def min_dem_rule(mod, i):
	return mod.dem[i]  >= mod.demmin[i]
mod.MinDemConstraint = Constraint(mod.N, rule=min_dem_rule)


#Unmet Demand is less than Demand
#	("Unmet Demand") <= demand
def unmet_rule(mod,i):
	return mod.unmet[i] <= mod.dem[i]
mod.UnmetConstraint = Constraint(mod.N, rule=unmet_rule)

#Unmet Demand Dual 	[Phi^(D.min), Phi^(D.max)]
def unmet_ruleZM(mod,i):
	return mod.dem[i] - mod.unmet[i] <= mod.M * mod.z[i]
mod.UnmetConstraintZM2 = Constraint(mod.N, rule=unmet_ruleZM)


#Transmisson Capacity
#	abs(transmission) <= capacity * lines
def cap_rule(mod, i, j):
	return mod.tran[i,j] <=  mod.cap[i,j] * mod.x_star[i,j]
mod.CapConstraint = Constraint(mod.L, rule=cap_rule)
def cap_rule2(mod, i, j):
	return mod.tran[i,j] >=  -mod.cap[i,j] * mod.x_star[i,j]
mod.CapConstraint2 = Constraint(mod.L, rule=cap_rule2)

#Transmission Capacicty Dual	[Phi^(L.Max) ; Phi^(L.Min)]
#	<= z*M
def cap_ruleZM1(mod, i, j):
	zmem = mod.LLen
	n = mod.NLen
	return (- mod.tran[i,j] + mod.cap[i,j] * mod.x_star[i,j] <=
		mod.M * mod.z[zmem - (n-i)*(n-i-1)/2 + j - i - 1])
mod.CapConstraintZM1 = Constraint(mod.L, rule=cap_ruleZM1)
def cap_ruleZM2(mod, i, j):
	zmem = 2*mod.LLen
	n = mod.NLen
	return (mod.tran[i,j] + mod.cap[i,j] * mod.x_star[i,j] <=
		mod.M *  mod.z[zmem - (n-i)*(n-i-1)/2 + j - i - 1])
mod.CapConstraintZM2= Constraint(mod.L, rule=cap_ruleZM2)


#Flow (Supply and Demand)
#Alpha is unfilled demand
#	Flow + Supply - Demand >= Unmet
#flow is positive if from Node 1 to 2, and negative if Node 2 to 1
def flow_rule(mod, i):
	flowcol = sum(mod.tran[j,j2] for (j,j2) in mod.L
		if (j2 == i)) 
	flowrow = sum(mod.tran[j,j2] for (j,j2) in mod.L
		if (j == i))
	return (flowcol - flowrow + mod.gen[i] - mod.dem[i]
		>=  -mod.unmet[i])
mod.FlowConstraint = Constraint(mod.N, rule=flow_rule)

#Flow Dual		[Lambda]
#	<=ZM
def flow_ruleZM(mod, i):
	zmem = 4*mod.LLen + 2*mod.NLen
	flowcol = sum(mod.tran[j,j2] for (j,j2) in mod.L
		if (j2 == i)) 
	flowrow = sum(mod.tran[j,j2] for (j,j2) in mod.L
		if (j == i))
	return (flowcol + mod.gen[i] - flowrow - mod.dem[i] +  mod.unmet[i]
		<=  mod.M * mod.z[zmem + i])
mod.FlowConstraintZM = Constraint(mod.N, rule=flow_ruleZM)


#Adjust the supply uncertainty budget
#	If supmin = supmax, no need to calculate, it is feasible
def unc_sup_rule(mod):
	if sum(mod.supmax[i] - mod.supmin[i] for i in mod.N) <= RC.TOL:
		return Constraint.Feasible
	else:
		return  (sum(mod.supmax[i] - mod.genpos[i] for i in mod.N)
			/ sum(mod.supmax[i] - mod.supmin[i] for i in mod.N)
			>= mod.uncS)
mod.UncSupConstraint = Constraint(rule=unc_sup_rule)


#Adjust the demand uncertainty budget
#	If demmin = demmax, no need to calculate, it is feasible
def unc_dem_rule(mod):
	if sum(mod.demmax[i] - mod.demmin[i] for i in mod.N) <= RC.TOL:
		return Constraint.Feasible
	else:
		return  (sum(mod.dem[i] - mod.demmin[i] for i in mod.N)
			/ sum(mod.demmax[i] - mod.demmin[i] for i in mod.N)
			<= mod.uncD)
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


# Theta Rules
# Theta is the angle at each node
# For each line from i to j:
#	|Transmision - (b)*(theta_i - theta_j)| <= M(1-route_on)
def theta_rule1(mod,i,j):
	return (mod.tran[i,j] - mod.b[i,j] * (mod.theta[i] - mod.theta[j])
		<= (1 - mod.route_on[i,j]) * mod.Mtheta)
mod.ThetaConstraint1 = Constraint(mod.L, rule=theta_rule1)
def theta_rule2(mod,i,j):
	return (-(mod.tran[i,j] - mod.b[i,j] * (mod.theta[i] - mod.theta[j]))
		<= (1 - mod.route_on[i,j]) * mod.Mtheta)
mod.ThetaConstraint2 = Constraint(mod.L, rule=theta_rule2)

# Theta Rules Dual		[Phi^L]
#
#

# -Pi <= Phi <= Phi Dual
#
#

'''
#KEEP?
#
def ref_rule1(mod,i):
	return (mod.theta[i] <= mod.ref[i]) 
mod.RefConstraint1 = Constraint(mod.N, rule=ref_rule1)
def ref_rule2(mod,i):
	return (-mod.theta[i] <= mod.ref[i]) 
mod.RefConstraint2 = Constraint(mod.N, rule=ref_rule2)
#
'''

 

##########
#TO TEST
##########


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


