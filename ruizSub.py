# -*- coding: utf-8 -*-
from pyomo.environ import *
from pyomo.opt import SolverFactory
import ruizC as RC


############################################################
#SubProblem

#Create an abstract model in Pyomo for the Sub Problem
#Input is expected to be .dat file
#Will take x_star (binary) as a parameter
############################################################

sub = AbstractModel()
opt = SolverFactory(RC.SOLVER)

#Parameters
sub.N = 	Set()						#Nodes
sub.L = 	Set(within=sub.N*sub.N)		#Lines
sub.c = 	Param(sub.L)				#cost per line
sub.cap =	Param(sub.L)				#line capacity
sub.sigma = Param()						#Hours in a year
sub.demmax = Param(sub.N)				#Maximum possible Demand
sub.demmin = Param(sub.N)				#Minimum possible Demand
sub.supmax = Param(sub.N)				#Maximum possible Supply
sub.supmin = Param(sub.N)				#Minimum possible Supply
sub.M	=	Param()						#Max M for Fourtuny-Amat
sub.gencost =	Param(sub.N)			#Cost to generate	
sub.shed =  Param(sub.N)				#Load Shedding Cost Per Node
sub.uncD =  Param()						#Uncertainty in Demand	
sub.uncS =  Param()						#Uncertainty in Supply
sub.x_star = Param(sub.L, domain=NonNegativeIntegers, default=0,
	mutable = True) 					#Built Lines

sub.conLen = Param()					#Constraints in Primal 
sub.constraints = RangeSet(1,sub.conLen) #(1, '# of constraints')
sub.varLen = Param()					#Variables in Primal
sub.NLen = Param()						#Length of N
sub.LLen = Param()						#Length of L
sub.LRange = RangeSet(1,sub.LLen)		#(1, '# of lines')

#B^T Matrix
#Split into parts for each variable
sub.gen_mu = Param(sub.N*sub.constraints)
sub.alpha_mu = Param(sub.N*sub.constraints)
sub.tran_mu = Param(RangeSet(1,sub.LLen)*sub.constraints)

#Variables
sub.tran = Var(sub.L, within=Reals) 			#Ammount Transmitted
sub.gen  =	Var(sub.N, domain=NonNegativeReals)	#Generation
sub.genpos = Var(sub.N, domain=NonNegativeReals)#Max Possible Generation
sub.alpha = Var(sub.N, domain=NonNegativeReals)	#Unfilled Demand
sub.dem  =	Var(sub.N, domain=NonNegativeReals)	#Demand
sub.lamda = Var(sub.constraints)				#Lamda Free
sub.mu	=	Var(sub.constraints, domain=NonNegativeReals)	#Mu >= 0
sub.z	=	Var(sub.constraints, domain=Binary)		#For Slackness

#############################################
#Functions
#############################################

#Objective Function
#	max [gen_cost * generation + shed_cost * unfilled_demand]
def obj_expression(sub):
	return (sub.sigma * (sum(sub.gencost[i] * sub.gen[i] for i in sub.N) 
		 + sum(sub.shed[i] * sub.alpha[i] for i in sub.N))
		 + sum(sub.c[j] * sub.x_star[j] for j in sub.L))
sub.OBJ = Objective(rule=obj_expression, sense = maximize)


#Transmisson Capacity
#	abs(transmission) <= capacity * lines
def cap_rule(sub, i, j):
	return sub.tran[i,j] <=  sub.cap[i,j] * sub.x_star[i,j]
sub.CapConstraint = Constraint(sub.L, rule=cap_rule)
def cap_rule2(sub, i, j):
	return sub.tran[i,j] >=  -sub.cap[i,j] * sub.x_star[i,j]
sub.CapConstraint2 = Constraint(sub.L, rule=cap_rule2)
#	<= z*M
def cap_rule_zm1(sub, i, j):
	zmem = sub.LLen
	n = sub.NLen
	return (- sub.tran[i,j] + sub.cap[i,j] * sub.x_star[i,j] <=
		sub.M * sub.z[zmem - (n-i)*(n-i-1)/2 + j - i - 1])
sub.CapConstraintZM1 = Constraint(sub.L, rule=cap_rule_zm1)
def cap_rule_zm2(sub, i, j):
	zmem = 2*sub.LLen
	n = sub.NLen
	return (sub.tran[i,j] + sub.cap[i,j] * sub.x_star[i,j] <=
		sub.M *  sub.z[zmem - (n-i)*(n-i-1)/2 + j - i - 1])
sub.CapConstraintZM2= Constraint(sub.L, rule=cap_rule_zm2)



#Unmet Demand is less than Demand
#	alpha("unmet demand") <= demand
def alpha_rule(sub,i):
	return sub.alpha[i] <= sub.dem[i]
sub.AlphaConstraint = Constraint(sub.N, rule=alpha_rule)
#	<= ZM
def alpha_ruleZM(sub,i):
	zmem = 4*sub.LLen + sub.NLen
	return sub.dem[i] - sub.alpha[i] <= sub.M * sub.z[zmem + i]
sub.AlphaConstraintZM2 = Constraint(sub.N, rule=alpha_ruleZM)


#Flow (Supply and Demand)
#Alpha is unfilled demand
#	flow + supply - demand >= alpha
#flow is positive if from Node 1 to 2, and negative if Node 2 to 1
#	alpha <= demand	
def flow_rule(sub, i):
	flowcol = sum(sub.tran[j,j2] for (j,j2) in sub.L
		if (j2 == i)) 
	flowrow = sum(sub.tran[j,j2] for (j,j2) in sub.L
		if (j == i))
	return (flowcol + sub.gen[i] - flowrow - sub.dem[i]
		>=  -sub.alpha[i])
sub.FlowConstraint = Constraint(sub.N, rule=flow_rule)
#	<=ZM
def flow_ruleZM(sub, i):
	zmem = 4*sub.LLen + 2*sub.NLen
	flowcol = sum(sub.tran[j,j2] for (j,j2) in sub.L
		if (j2 == i)) 
	flowrow = sum(sub.tran[j,j2] for (j,j2) in sub.L
		if (j == i))
	return (flowcol + sub.gen[i] - flowrow - sub.dem[i] +  sub.alpha[i]
		<=  sub.M * sub.z[zmem + i])
sub.FlowConstraintZM = Constraint(sub.N, rule=flow_ruleZM)

#Supply min and max
#	generation <= possible_gen <= gen_max
def genpos_rule(sub, i):
	return (sub.genpos[i] <= sub.supmax[i])
sub.GenPosConstraint = Constraint(sub.N, rule=genpos_rule)
def max_sup_rule(sub, i):
	return (sub.gen[i] <= sub.genpos[i])
sub.MaxSupConstraint = Constraint(sub.N, rule=max_sup_rule)

#def min_sup_rule(sub, i):
#	return (sub.supmin[i] <= sub.gen[i])
#sub.MinSupConstraint = Constraint(sub.N, rule=min_sup_rule)


#Demand Min and Max
#	demand_min <= demand <= demand_max
def max_dem_rule(sub, i):
	return sub.dem[i]  <= sub.demmax[i]
sub.MaxDemConstraint = Constraint(sub.N, rule=max_dem_rule)
def min_dem_rule(sub, i):
	return sub.dem[i]  >= sub.demmin[i]
sub.MinDemConstraint = Constraint(sub.N, rule=min_dem_rule)


#Adjust the supply uncertainty budget
def unc_sup_rule(sub):
	if sum(sub.supmax[i] - sub.supmin[i] for i in sub.N) <= RC.TOL:
		return Constraint.Feasible
	else:
		return  (sum(sub.supmax[i] - sub.genpos[i] for i in sub.N)
			/ sum(sub.supmax[i] - sub.supmin[i] for i in sub.N)
			== sub.uncS)
sub.UncSupConstraint = Constraint(rule=unc_sup_rule)


#Adjust the demand uncertainty budget
def unc_dem_rule(sub):
	if sum(sub.demmax[i] - sub.demmin[i] for i in sub.N) <= RC.TOL:
		return Constraint.Feasible
	else:
		return  (sum(sub.dem[i] - sub.demmin[i] for i in sub.N)
			/ sum(sub.demmax[i] - sub.demmin[i] for i in sub.N)
			>= sub.uncD)
sub.UncDemConstraint = Constraint(rule=unc_dem_rule)


#Lambda Mu Yp Rules
#############################################
#Generation
#I = 6
def lamda_mu_gen_rule(sub, i):
	return (sum(sub.gen_mu[i,j] * sub.mu[j] for j in sub.constraints)
		== -sub.gencost[i])
sub.LamdaMuGenConstraint = Constraint(sub.N, rule=lamda_mu_gen_rule)

#Alpha
# I = 6
def lamda_mu_alpha_rule(sub, i):
	return (sum(sub.alpha_mu[i,j] * sub.mu[j] for j in sub.constraints)
		== -sub.shed[i])
sub.LamdaMuAlphaConstraint = Constraint(sub.N, rule=lamda_mu_alpha_rule)

#Transmission
# I = 15
def lamda_mu_tran_rule(sub, i):
	return (sum(sub.tran_mu[i,j] * sub.mu[j] for j in sub.constraints)
		== 0)
sub.LamdaMuTranConstraint = Constraint(sub.LRange, rule=lamda_mu_tran_rule)
#############################################


#Mu is 0 when z is 1	
def mu_bigm_rule(sub, i):
	return sub.mu[i] <= sub.M*(1-sub.z[i]) 
sub.MuBigM = Constraint(sub.constraints, rule=mu_bigm_rule)

 

##########
#TO TEST
##########
'''
isub = sub.create_instance(RC.SUB)
results = opt.solve(isub, tee=True)
#results = opt.solve(isub)
#isub.pprint()
results.write()

##To Print

for v in isub.component_objects(Var, active=True):
	print ("Variable",v)
	varob = getattr(isub, str(v))
	for index in varob:
		print ("   ",index, varob[index].value) 
'''

