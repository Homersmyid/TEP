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
opt.options['mipgap'] = RC.MIPGAP

###############################################################
#Parameters and Sets
###############################################################

#Sets
mod.N = 	Set()						#Nodes
mod.L = 	Set(within=mod.N*mod.N)		#Lines

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

#Parameters that come from Master
mod.x_star = 	Param(mod.L, domain=NonNegativeIntegers, default=0,
	mutable = True) 					#Built Lines


###############################################################
#Variables
###############################################################

#Variables
mod.tran   = Var(mod.L, within=Reals) 			 #Ammount Transmitted
mod.dem    = Var(mod.N, domain=NonNegativeReals) #Demand
mod.genpos = Var(mod.N, domain=NonNegativeReals) #Max Possible Gen
mod.gen    = Var(mod.N, domain=NonNegativeReals) #Generation
mod.unmet  = Var(mod.N, domain=NonNegativeReals) #Unfilled Demand

#Angle in [-pi,pi]
mod.theta = Var(mod.N,bounds=(-math.pi, math.pi))

#Dual Variables
#Positive for Inequalities
mod.genmax_dual	  =	Var(mod.N, domain=NonNegativeReals) #[Phi^(E.max)]
mod.genmin_dual	  =	Var(mod.N, domain=NonNegativeReals) #[Phi^(E.min)]
mod.unmetmax_dual = Var(mod.N, domain=NonNegativeReals) #[Phi^(D.max)]
mod.unmetmin_dual = Var(mod.N, domain=NonNegativeReals) #[Phi^(D.min)]
mod.thetamax_dual = Var(mod.N, domain=NonNegativeReals) #[Phi^(N.max)]
mod.thetamin_dual = Var(mod.N, domain=NonNegativeReals) #[Phi^(N.min)]
mod.capmax_dual   =	Var(mod.L, domain=NonNegativeReals) #[Phi^(L.max)]
mod.capmin_dual   =	Var(mod.L, domain=NonNegativeReals) #[Phi^(L.min)]
#Free for Equalities (thus no complementarity)
mod.theta_dual 	  = Var(mod.L) 							#[Phi^(L)]
mod.flow_dual 	  =	Var(mod.N) 							#[Lambda]
mod.ref_dual	  = Var()								#[Chi^Ref]
	#Dual for reference variable -->  theta[mod.ref] = 0

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
#	max sigma * [gen_cost * generation + shed_cost * unfilled_demand]
#		+ c^t*x
def obj_expression(mod):
	return (mod.sigma * (sum(mod.gencost[i] * mod.gen[i] for i in mod.N) 
		 + sum(mod.shed[i] * mod.unmet[i] for i in mod.N))
		 + sum(mod.c[j] * mod.x_star[j] for j in mod.L))
mod.Obj = Objective(rule=obj_expression, sense = maximize)


#Gen Possible Min and Max
# Gen_pos is generation possible at a node
#	gen_min <= possible_gen <= gen_max
def genpos_rule(mod, i):
	return (mod.supmin[i], mod.genpos[i], mod.supmax[i])
mod.GenPosConstraint = Constraint(mod.N, rule=genpos_rule)

	
#Generation is less than Gen Possible
# Genation is the decision variable being minimized by the agent
# Generation <= Generation Possible
def max_gen_rule(mod, i):
	return (mod.gen[i] <= mod.genpos[i])
mod.MaxGenConstraint = Constraint(mod.N, rule=max_gen_rule)


#Demand Min and Max
#	demand_min <= (Demand) <= demand_max
def max_dem_rule(mod, i):
	return (mod.demmin[i], mod.dem[i], mod.demmax[i])
mod.MaxDemConstraint = Constraint(mod.N, rule=max_dem_rule)


#Unmet Demand is less than Demand
#	("Unmet Demand") <= demand
def unmet_rule(mod,i):
	return mod.unmet[i] <= mod.dem[i]
mod.UnmetConstraint = Constraint(mod.N, rule=unmet_rule)


#Transmisson Capacity
#	abs(transmission) <= capacity * lines
def cap_rule(mod, i, j):
	return (-mod.cap[i,j] * mod.x_star[i,j], mod.tran[i,j],
		mod.cap[i,j] * mod.x_star[i,j]) 
mod.CapConstraint = Constraint(mod.L, rule=cap_rule)


#Flow (Supply and Demand)
#Alpha is unfilled demand
#	Flow + Supply - Demand = (-Unmet)
#flow is positive if from Node 1 to 2, and negative if Node 2 to 1
def flow_rule(mod, i):
	flowcol = sum(mod.tran[j,j2] for (j,j2) in mod.L if (j2 == i))
	flowrow = sum(mod.tran[j,j2] for (j,j2) in mod.L if (j == i))
	return flowcol - flowrow + mod.gen[i] - mod.dem[i] == -mod.unmet[i]
mod.FlowConstraint = Constraint(mod.N, rule=flow_rule)


# Theta Rules
# Theta is the angle at each node
# For each line from i to j:
#	(b)*(theta_i - theta_j) = flow  (if route "ij" active)
def theta_rule(mod,i,j):
	return ( (mod.b[i,j] * (mod.theta[i] - mod.theta[j]))
		* mod.x_star[i,j] == mod.tran[i,j])
mod.ThetaConstraint = Constraint(mod.L, rule=theta_rule)



# Reference Theta
#	Theta_refernce = 0
def ref_rule(mod):
	return mod.theta[value(mod.ref)] == 0
mod.RefConstraint = Constraint(rule=ref_rule)

###############################
#Uncertainty Budgets
###############################

#Adjust the supply uncertainty budget
#	If supmin = supmax, no need to calculate, it is feasible
def unc_sup_rule(mod):
	if (sum(mod.supmax[i] - mod.supmin[i] for i in mod.N)
		<= (RC.UNCTOL * len(mod.N))):
		return Constraint.Feasible
	else:
		return  (sum(mod.supmax[i] - mod.genpos[i] for i in mod.N)
			/ sum(mod.supmax[i] - mod.supmin[i] for i in mod.N)
			== mod.uncS)
mod.UncSupConstraint = Constraint(rule=unc_sup_rule)


#Adjust the demand uncertainty budget
#	If demmin = demmax, no need to calculate, it is feasible
def unc_dem_rule(mod):
	if (sum(mod.demmax[i] - mod.demmin[i] for i in mod.N)
		<= (RC.UNCTOL * len(mod.N))):
		return Constraint.Feasible
	else:
		return  (sum(mod.dem[i] - mod.demmin[i] for i in mod.N)
			/ sum(mod.demmax[i] - mod.demmin[i] for i in mod.N)
			== mod.uncD)
mod.UncDemConstraint = Constraint(rule=unc_dem_rule)

#############################################
#Linearized Complementarity Constraints
#############################################

#Generation Max Dual	[Phi^(E.max)]
def gen_rule_max_dual1(mod,i):
	return  mod.genpos[i] - mod.gen[i] <= (mod.Mgen * mod.z_genmax[i]) 
mod.GenMaxConstraintDual1 = Constraint(mod.N, rule=gen_rule_max_dual1)
def gen_rule_max_dual2(mod,i):
	return mod.genmax_dual[i] <= (mod.M * (1 - mod.z_genmax[i])) 
mod.GenMaxConstraintDual2 = Constraint(mod.N, rule=gen_rule_max_dual2)

#Generation Min Dual	[Phi^(E.min)]
def gen_rule_min_dual1(mod,i):
	return mod.gen[i] <= mod.Mgen * mod.z_genmin[i] 
mod.GenMinConstraintDual1 = Constraint(mod.N, rule=gen_rule_min_dual1)
def gen_rule_min_dual2(mod,i):
	return mod.genmin_dual[i] <= (mod.M * (1 - mod.z_genmin[i])) 
mod.GenMinConstraintDual2 = Constraint(mod.N, rule=gen_rule_min_dual2)

###########

#Unmet Demand Max Dual 	[Phi^(D.max)]
def unmet_rule_max_dual1(mod,i):
	return mod.dem[i] - mod.unmet[i] <= mod.Mdem * mod.z_unmetmax[i] 
mod.UnmetMaxConstraint1 = Constraint(mod.N, rule=unmet_rule_max_dual1)
def unmet_rule_max_dual2(mod,i):
	return mod.unmetmax_dual[i] <= mod.M * (1 - mod.z_unmetmax[i]) 
mod.UnmetMaxConstraint2 = Constraint(mod.N, rule=unmet_rule_max_dual2)

#Unmet Demand Min Dual 	[Phi^(D.min)]
def unmet_rule_min_dual1(mod,i):
	return mod.unmet[i] <= mod.Mdem * mod.z_unmetmin[i] 
mod.UnmetMinConstraint1 = Constraint(mod.N, rule=unmet_rule_min_dual1)
def unmet_rule_min_dual2(mod,i):
	return mod.unmetmin_dual[i] <= mod.M * (1 - mod.z_unmetmin[i]) 
mod.UnmetMinConstraint2 = Constraint(mod.N, rule=unmet_rule_min_dual2)

###########

#Transmission Capacicty Max Dual	[Phi^(L.Max)]
def cap_rule_max_dual1(mod, i, j):
	return ( (mod.cap[i,j] * mod.x_star[i,j]) - mod.tran[i,j]
		<= (mod.Mcap * mod.z_capmax[i,j]) )
mod.CapMaxConstraintDual1 = Constraint(mod.L, rule=cap_rule_max_dual1)
def cap_rule_max_dual2(mod, i, j):
	return mod.capmax_dual[i,j] <= (mod.M * (1 - mod.z_capmax[i,j])) 
mod.CapMaxConstraintDual2 = Constraint(mod.L, rule=cap_rule_max_dual2)

#Transmission Capacicty Min Dual	[Phi^(L.Min)]
def cap_rule_min_dual1(mod, i, j):
	return ( (mod.cap[i,j] * mod.x_star[i,j]) + mod.tran[i,j]
		<= mod.Mcap * mod.z_capmin[i,j] )
mod.CapMinConstraintDual1 = Constraint(mod.L, rule=cap_rule_min_dual1)
def cap_rule_min_dual2(mod, i, j):
	return mod.capmin_dual[i,j] <= (mod.M * (1 - mod.z_capmin[i,j])) 
mod.CapMinConstraintDual2 = Constraint(mod.L, rule=cap_rule_min_dual2)

###########

# Theta Rules Max Dual		[Phi^(N.max)]
def theta_rule_max_dual1(mod,i):
	return  math.pi - mod.theta[i] <= mod.Mtheta * mod.z_thetamax[i] 
mod.ThetaMaxConstraint1 = Constraint(mod.N, rule=theta_rule_max_dual1)
def theta_rule_max_dual2(mod,i):
	return mod.thetamax_dual[i] <= mod.M * (1 - mod.z_thetamax[i]) 
mod.ThetaMaxConstraint2 = Constraint(mod.N, rule=theta_rule_max_dual2)

# Theta Rules Min Dual		[Phi^(N.min)]
def theta_rule_min_dual1(mod,i):
	return mod.theta[i] + math.pi <= mod.Mtheta * mod.z_thetamin[i] 
mod.ThetaMinConstraint1 = Constraint(mod.N, rule=theta_rule_min_dual1)
def theta_rule_min_dual2(mod,i):
	return mod.thetamin_dual[i] <= mod.M * (1 - mod.z_thetamin[i]) 
mod.ThetaMinConstraint2 = Constraint(mod.N, rule=theta_rule_min_dual2)

##################################
#Differentiating the Lagrangian
##################################

#Lagrangian with respect to Generation
#	(Sigma * GenCost) - Flow_Dual + GenMax_Dual - GenMin_Dual = 0
def lag_gen(mod,i):
	return ((mod.sigma * mod.gencost[i]) - mod.flow_dual[i] 
			+ mod.genmax_dual[i] -  mod.genmin_dual[i] == 0)
mod.LagrangianGenConstraint = Constraint(mod.N, rule=lag_gen)


#Lagrangian with respect to Unmet Demand
#	(Sigma * Unmet_Penalty) - Flow_dual + UnmetMax_Dual -UnmetMin_Dual=0
def lag_unmet(mod,i):
	return ((mod.sigma * mod.shed[i]) - mod.flow_dual[i]
			+ mod.unmetmax_dual[i] -  mod.unmetmin_dual[i] == 0) 
mod.LagrangianUnmetConstraint = Constraint(mod.N, rule=lag_unmet)


#Lagrangian with respect to Tranmission
#	For each line in network from source --> dest:
#		Flow_Dual(source) - Flow_Dual(dest) - Theta_Dual
#		+ Capmax_Dual - Capmin_Dual = 0
def lag_trans(mod,i,j):
	return ( (mod.flow_dual[i] - mod.flow_dual[j] - mod.theta_dual[i,j]
				+ mod.capmax_dual[i,j] - mod.capmin_dual[i,j] )
			 * mod.x_star[i,j] == 0)
mod.LagrangianTransConstraint = Constraint(mod.L, rule=lag_trans)


#Lagrangian with respect to Theta
#	At each node:
#		Sum(B*x*Theta_Dual)_(Send) -  Sum(B*x*Theta_Dual)_(Recieve)
#		+ ThetaMax_Dual + ThetaMin_Dual = 0
#Only consider Lines that are in network
def lag_theta(mod,i):	
	if_ref = 1 if i == mod.ref else 0
	return (  sum(mod.x_star[j,j2] * mod.b[j,j2] * mod.theta_dual[j,j2]
					for (j,j2) in mod.L if (j == i))
		    - sum(mod.x_star[j,j2] * mod.b[j,j2] * mod.theta_dual[j,j2]
					for (j,j2) in mod.L if (j2 == i))
			+ mod.thetamax_dual[i] - mod.thetamin_dual[i]
			- (mod.ref_dual * if_ref) == 0 )
mod.LagrangianThetaConstraint = Constraint(mod.N, rule=lag_theta)	


 

##########
#TO TEST
##########
'''
isub = mod.create_instance(RC.DATA)

#Lines At Start
START_X_STAR = [(1,2,1), (1,3,1), (1,4,0), (1,5,0), (1,6,0),
				(2,3,0), (2,4,0), (2,5,0), (2,6,0), (3,4,0),
				(3,5,0), (3,6,0), (4,5,1), (4,6,0), (5,6,0)]

#Set x_star 
for xi in START_X_STAR:
	isub.x_star[xi[0], xi[1]] = xi[2]

results = opt.solve(isub)
results.write()


##To Print
for v in isub.component_objects(Var, active=True):
	print ("Variable",v)
	varob = getattr(isub, str(v))
	for index in varob:
		print ("   ",index, varob[index].value) 
#isub.pprint()


'''



'''
file1 = open("junk", 'w') 
results.write()
##To Print
for v in isub.component_objects(Var, active=True):
	file1.write("\n")
	file1.write("Variable ")
	file1.write(str(v))
	file1.write("\n")
	varob = getattr(isub, str(v))
	for index in varob:
		file1.write("   ")
		file1.write(str(index))
		file1.write("   ")
		file1.write(str(varob[index].value)) 
		file1.write("\n")  
isub.pprint()
file1.close()
'''
