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

#Parameters that come from Master
mod.x_star = 	Param(mod.L, domain=NonNegativeIntegers, default=0,
	mutable = True) 					#Built Lines


###############################################################
#Variables
###############################################################

#Variables
mod.tran = 	Var(mod.L, within=Reals) 			#Ammount Transmitted
mod.gen  =	Var(mod.N, domain=NonNegativeReals)	#Generation
mod.genpos = Var(mod.N, domain=NonNegativeReals)#Max Possible Generation
mod.unmet = Var(mod.N, domain=NonNegativeReals)	#Unfilled Demand
mod.dem  =	Var(mod.N, domain=NonNegativeReals)	#Demand

#Angle in [-pi,pi]
mod.theta = Var(mod.N,bounds=(-math.pi, math.pi), initialize=0)

#If route used (Bianary). Needed for when the lines per route is > 1
mod.route_on =	Var(mod.L, domain=Binary, initialize=0)		

#Dual Variables
mod.genmax_dual	=	Var(mod.N, domain=NonNegativeReals) #[Phi^(E.max)]
mod.genmin_dual	=	Var(mod.N, domain=NonNegativeReals) #[Phi^(E.min)]
mod.unmetmax_dual = Var(mod.N, domain=NonNegativeReals) #[Phi^(D.max)]
mod.unmetmin_dual = Var(mod.N, domain=NonNegativeReals) #[Phi^(D.min)]
mod.thetamax_dual = Var(mod.N, domain=NonNegativeReals) #[Phi^(N.max)]
mod.thetamin_dual = Var(mod.N, domain=NonNegativeReals) #[Phi^(N.min)]
mod.capmax_dual = 	Var(mod.L, domain=NonNegativeReals) #[Phi^(L.max)]
mod.capmin_dual = 	Var(mod.L, domain=NonNegativeReals) #[Phi^(L.min)]
mod.flow_dual	=	Var(mod.L, domain=NonNegativeReals) #[Lambda]
	#Lamda is the only dual for an equality, thus no complementarity


#Binary Variables for Complementary Condition Linearization
mod.z_genmax	=	Var(mod.N, domain=Binary)	#For Generation Max
mod.z_genmin	=	Var(mod.N, domain=Binary)	#For Generation Min
mod.z_unmetmax	=	Var(mod.N, domain=Binary)	#For Unmet Demand Max
mod.z_unmetmin	=	Var(mod.N, domain=Binary)	#For Unmet Demand Min
mod.z_thetamax	=	Var(mod.N, domain=Binary)	#For Theta Max 
mod.z_thetamin	=	Var(mod.N, domain=Binary)	#For Theta Min
mod.z_capmax	=	Var(mod.L, domain=Binary)	#For Capacity Max
mod.z_capmin	=	Var(mod.L, domain=Binary)	#For Capacity Min


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


# Route Rule
# To see what routes are activated, route_on is binary
# MaxLines * Route_on >= x
def route_rule(mod, i,j):
	return (mod.maxLines * mod.route_on[i,j] >= mod.x_star[i,j])
mod.RouteConstraint = Constraint(mod.L, rule=route_rule)

##############################################################

#GenPos Min and Max
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
mod.MaxGenConstraint = Constraint(mod.N, rule=max_gen_rule)

#Generation Max Dual	[Phi^(E.max)]
def gen_rule_max_dual1(mod,i):
	return  mod.genpos[i] - mod.gen[i] <= (mod.M * mod.z_genmax[i])
mod.GenMaxConstraintDual1 = Constraint(mod.N, rule=gen_rule_max_dual1)
def gen_rule_max_dual2(mod,i):
	return mod.genmax_dual[i] <= (mod.M * (1 - mod.z_genmax[i]))
mod.GenMaxConstraintDual2 = Constraint(mod.N, rule=gen_rule_max_dual2)

#Generation Min Dual	[Phi^(E.min)]
def gen_rule_min_dual1(mod,i):
	return mod.gen[i] <= mod.M * mod.z_genmin[i]
mod.GenMinConstraintDual1 = Constraint(mod.N, rule=gen_rule_min_dual1)
def gen_rule_min_dual2(mod,i):
	return mod.genmin_dual[i] <= (mod.M * (1 - mod.z_genmin[i]))
mod.GenMinConstraintDual2 = Constraint(mod.N, rule=gen_rule_min_dual2)

##############################################################

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


#Unmet Demand Max Dual 	[Phi^(D.max)]
def unmet_rule_max_dual1(mod,i):
	return mod.dem[i] - mod.unmet[i] <= mod.M * mod.z_unmetmax[i]
mod.UnmetMaxConstraint1 = Constraint(mod.N, rule=unmet_rule_max_dual1)
def unmet_rule_max_dual2(mod,i):
	return mod.unmetmax_dual[i] <= mod.M * (1 - mod.z_unmetmax[i])
mod.UnmetMaxConstraint2 = Constraint(mod.N, rule=unmet_rule_max_dual2)


#Unmet Demand Min Dual 	[Phi^(D.min)]
def unmet_rule_min_dual1(mod,i):
	return mod.unmet[i] <= mod.M * mod.z_unmetmin[i]
mod.UnmetMinConstraint1 = Constraint(mod.N, rule=unmet_rule_min_dual1)
def unmet_rule_min_dual2(mod,i):
	return mod.unmetmin_dual[i] <= mod.M * (1 - mod.z_unmetmin[i])
mod.UnmetMinConstraint2 = Constraint(mod.N, rule=unmet_rule_min_dual2)

##############################################################

#Transmisson Capacity
#	abs(transmission) <= capacity * lines
def cap_rule(mod, i, j):
	return mod.tran[i,j] <=  mod.cap[i,j] * mod.x_star[i,j]
mod.CapConstraint = Constraint(mod.L, rule=cap_rule)
def cap_rule2(mod, i, j):
	return mod.tran[i,j] >=  -mod.cap[i,j] * mod.x_star[i,j]
mod.CapConstraint2 = Constraint(mod.L, rule=cap_rule2)


#Transmission Capacicty Max Dual	[Phi^(L.Max)]
def cap_rule_max_dual1(mod, i, j):
	return (mod.cap[i,j] * mod.x_star[i,j] - mod.tran[i,j]
		<= mod.M * mod.z_capmax[i,j])
mod.CapMaxConstraintDual1 = Constraint(mod.L, rule=cap_rule_max_dual1)
def cap_rule_max_dual2(mod, i, j):
	return mod.capmax_dual[i,j] <= (mod.M * (1 - mod.z_capmax[i,j]))
mod.CapMaxConstraintDual2 = Constraint(mod.L, rule=cap_rule_max_dual2)

#Transmission Capacicty Min Dual	[Phi^(L.Min)]
def cap_rule_min_dual1(mod, i, j):
	return (mod.cap[i,j] * mod.x_star[i,j] + mod.tran[i,j]
		<= mod.M * mod.z_capmin[i,j])
mod.CapMinConstraintDual1 = Constraint(mod.L, rule=cap_rule_min_dual1)
def cap_rule_min_dual2(mod, i, j):
	return mod.capmin_dual[i,j] <= (mod.M * (1 - mod.z_capmin[i,j]))
mod.CapMinConstraintDual2 = Constraint(mod.L, rule=cap_rule_min_dual2)

##############################################################

#Flow (Supply and Demand)
#Alpha is unfilled demand
#	Flow + Supply - Demand >= Unmet
#flow is positive if from Node 1 to 2, and negative if Node 2 to 1
def flow_rule(mod, i):
	flowcol = sum(mod.tran[j,j2] for (j,j2) in mod.L if (j2 == i))
	flowrow = sum(mod.tran[j,j2] for (j,j2) in mod.L if (j == i))
	return (flowcol - flowrow + mod.gen[i] - mod.dem[i]
		>=  -mod.unmet[i])
mod.FlowConstraint = Constraint(mod.N, rule=flow_rule)


##############################################################

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

##############################################################

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

# Theta Rules Max Dual		[Phi^(N.max)]
def theta_rule_max_dual1(mod,i):
	return  math.pi - mod.theta[i] <= mod.M * mod.z_thetamax[i]
mod.ThetaMaxConstraint1 = Constraint(mod.N, rule=theta_rule_max_dual1)
def theta_rule_max_dual2(mod,i):
	return mod.thetamax_dual[i] <= mod.M * (1 - mod.z_thetamax[i])
mod.ThetaMaxConstraint2 = Constraint(mod.N, rule=theta_rule_max_dual2)

# Theta Rules Min Dual		[Phi^(N.min)]
def theta_rule_min_dual1(mod,i):
	return mod.theta[i] + math.pi <= mod.M * mod.z_thetamin[i]
mod.ThetaMinConstraint1 = Constraint(mod.N, rule=theta_rule_min_dual1)
def theta_rule_min_dual2(mod,i):
	return mod.thetamin_dual[i] <= mod.M * (1 - mod.z_thetamin[i])
mod.ThetaMinConstraint2 = Constraint(mod.N, rule=theta_rule_min_dual2)

##############################################################

##################################
#Differentiating the Lagrangian
##################################


#Lagrangian with respect to Generation
#	(Sigma * GenCost) - Flow_Dual + GenMax_Dual - GenMin_Dual = 0
def lag_gen(mod,i):
	lambda_g = (sum(mod.flow_dual[j,j2] for (j,j2) in mod.L
		if ((j == i) or (j2 == i))))
	return ((mod.sigma * mod.gencost[i]) - lambda_g 
		+ mod.genmax_dual[i] -  mod.genmin_dual[i] == 0)
mod.LagrangianGenConstraint = Constraint(mod.N, rule=lag_gen)


#Lagrangian with respect to Unmet Demand
#	(Sigma * Unmet_Penalty) - Flow_dual + UnmetMax_Dual -UnmetMin_Dual=0
def lag_unmet(mod,i):
	lambda_d = sum(mod.flow_dual[j,j2] for (j,j2) in mod.L
		if ((j == i) or (j2 == i)))
	return ((mod.sigma * mod.shed[i]) - lambda_d 
		+ mod.unmetmax_dual[i] -  mod.unmetmin_dual[i] == 0)
mod.LagrangianUnmetConstraint = Constraint(mod.N, rule=lag_unmet)

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

#Set x_star 
for xi in RC.START_X_STAR:
	isub.x_star[xi[0], xi[1]] = xi[2]

results = opt.solve(isub, tee=True)
#isub.pprint()
results.write()

##To Print

for v in isub.component_objects(Var, active=True):
	print ("Variable",v)
	varob = getattr(isub, str(v))
	for index in varob:
		print ("   ",index, varob[index].value) 


