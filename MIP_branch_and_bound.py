# 0-1混合整数问题的分支定界算法
# wang

# creat LP
from gurobipy import *
import numpy as np
import copy
import matplotlib.pyplot as plt

RLP = Model('relaxed MIP')

#定义决策变量
rtsr = {}
rtmf = {}
rts = {}
rta = {}
qm = {}
qs = {}
ym = {}
ys = {}
dm = {}
dp = {}
for i in range(2):
    rtsr[i]={}
    for j in range(2):
        rtsr[i][j] = RLP.addVar(vtype=GRB.CONTINUOUS, name='rtsr_' + str(i)+ str(j))
for i in range(2):
    rtmf[i]={}
    for j in range(2):
        rtmf[i][j] = RLP.addVar(vtype=GRB.CONTINUOUS, name='rtmf_' + str(i)+ str(j))
for i in range(2):
    rts[i]={}
    for j in range(2):
        rts[i][j] = RLP.addVar(vtype=GRB.CONTINUOUS, name='rts_' + str(i)+ str(j))
for i in range(2):
    rta[i]={}
    for j in range(2):
        rta[i][j] = RLP.addVar(vtype=GRB.CONTINUOUS, name='rta_' + str(i)+ str(j))
for i in range(2):
    qm[i] = RLP.addVar(vtype=GRB.CONTINUOUS, name = 'qm_' + str(i))
for i in range(2):
    qs[i] = RLP.addVar(vtype=GRB.CONTINUOUS, name = 'qs_' + str(i))
for i in range(2):
    ym[i] = RLP.addVar(lb=0, ub=1, vtype = GRB.CONTINUOUS, name = 'ym_' + str(i))  
for i in range(2):
    ys[i] = RLP.addVar(lb=0, ub=1, vtype = GRB.CONTINUOUS, name = 'ys_' + str(i))
for i in range(2):
    dp[i] = RLP.addVar(vtype=GRB.CONTINUOUS, name='dp_' + str(i))
for i in range(2):
    dm[i] = RLP.addVar(vtype=GRB.CONTINUOUS, name='dm_' + str(i))

w = RLP.addVars(16,  vtype=GRB.CONTINUOUS, name='w')
z = RLP.addVars(16,  vtype=GRB.CONTINUOUS, name='z')
v = RLP.addVars(8,  vtype=GRB.CONTINUOUS, name='v')

#定义参数
wsr=800
wsf=900
wmf=800
dstsr = [ [459,585],[368,892] ]
dstsf = [ [900,1105],[1302,693] ]	
dstmf = [ [362,578],[605,759] ]
dmr = [500,544]
dmf = [1200,1265]

cs = [990,780]
cm = [1000,920]

fcs = [0.2,0.5]
fcm = [0.3,0.3]

eas = [13.86,10.92]
eam = [14, 12.88]
g = 300
k = 0.14
etr = 0.05
etra = 0.15
etrs = 0.03
rc = 1.46
tpa = 2.2
tps = 1.58
tb = 0.1

z1 = 1000000
z2 = 50000
P1 = 1000
P2 = 1

theta = 60

va = 0.08
vs = 1.2
T = 72

#设置目标
RLP.setObjective(-P1*dm[0]-P2*dp[1], GRB.MAXIMIZE)

#设置约束
RLP.addConstrs((qs[s] <= cs[s]*ys[s] for s in range(2)), name = 'ConA')
RLP.addConstrs((qm[m] <= cm[m]*ym[m] for m in range(2)), name = 'ConB')
RLP.addConstrs((qs[s] == gurobipy.quicksum(rtsr[s][r] for r in range(2)) + gurobipy.quicksum(rts[s][f]+rta[s][f] for f in range(2))for s in range(2)), name = 'ConC')
RLP.addConstrs((qm[m] == gurobipy.quicksum(rtmf[m][f] for f in range(2)) for m in range(2)), name ='ConD')
RLP.addConstrs((dmr[r] == gurobipy.quicksum(rtsr[s][r] for s in range(2)) for r in range(2)), name = 'ConE')
RLP.addConstrs((dmr[f] == gurobipy.quicksum((rts[s][f]+rta[s][f]) for s in range(2)) + gurobipy.quicksum(rtmf[m][f] for m in range(2)) for f in range(2)), name = 'ConF')
RLP.addConstr(gurobipy.quicksum(wsr*rtsr[s][r] for s in range(2) for r in range(2)) 
              - gurobipy.quicksum(ys[s]*fcs[s]-g*qs[s] for s in range(2)) - gurobipy.quicksum(rc*rtsr[s][r] for s in range(2) for r in range(2))
              + gurobipy.quicksum(wsf*(rts[s][f]+rta[s][f])-tps*rts[s][f]-tpa*rta[s][f]-tb*(rts[s][f]+rta[s][f]) for s in range(2) for f in range(2))
              + gurobipy.quicksum(wmf*rtmf[m][f] for m in range(2) for f in range(2))-gurobipy.quicksum(ym[m]*fcm[m]-g*qm[m]for m in range(2))
              - gurobipy.quicksum(rtmf[m][f]*rc for m in range(2) for f in range(2))
              +dm[0]-dp[0] == z1, name='ConG')
RLP.addConstr(gurobipy.quicksum(etr*dstsr[s][r]*rtsr[s][r] for s in range(2) for r in range(2))
              + gurobipy.quicksum(etr*dstmf[m][f]*rtmf[m][f]+dm[1]-dp[1] for m in range(2) for f in range(2))
              + gurobipy.quicksum(0.5*w[i] + z[i] for i in range(16)) + gurobipy.quicksum((v[j]*v[j])/(2*theta) for j in range(8))
              <= - gurobipy.quicksum(eas[s]*ys[s]+k*qs[s] for s in range(2)) - gurobipy.quicksum(eam[m]*ym[m]+k*qm[m] for m in range(2))
              - gurobipy.quicksum(etra*dstsf[s][f]*rta[s][f]+etrs*dstsf[s][f]*rts[s][f] for s in range(2) for f in range(2)) 
              + z2, name = 'ConH')
RLP.addConstrs((w[i]-w[i+8] == 0.005*v[i] for i in range(8)), name = 'ConI')
RLP.addConstr(z[0]-z[8] == 0.005*dstsr[0][0]*rtsr[0][0], name = 'ConJ')
RLP.addConstr(z[1]-z[9] == 0.005*dstsr[0][1]*rtsr[0][1], name = 'ConK')
RLP.addConstr(z[2]-z[10] == 0.005*dstsr[1][0]*rtsr[1][0], name = 'ConL')
RLP.addConstr(z[3]-z[11] == 0.005*dstsr[1][1]*rtsr[1][1], name = 'ConM')
RLP.addConstr(z[4]-z[12] == 0.005*dstmf[0][0]*rtmf[0][0], name = 'ConN')
RLP.addConstr(z[5]-z[13] == 0.005*dstmf[0][1]*rtmf[0][1], name = 'ConO')
RLP.addConstr(z[6]-z[14] == 0.005*dstmf[1][0]*rtmf[1][0], name = 'ConP')
RLP.addConstr(z[7]-z[15] == 0.005*dstmf[1][1]*rtmf[1][1], name = 'ConQ')
RLP.addConstr(gurobipy.quicksum(va*rta[s][f]+vs*rts[s][f] for s in range(2) for f in range(2))<=T, name='ConR')

RLP.optimize()


# Node class
class Node:
# this class defines the node
    def __init__(self):       #初始化节点的设置
        self.local_LB = -5000
        self.local_UB = 0
        self.x_sol = {}
        self.x_int_sol = {}
        self.branch_var_list = []
        self.model = None
        self.cnt = None
        self.is_integer = False

    def deepcopy_node(node):   #复制节点
        new_node = Node()
        new_node.local_LB = -5000
        new_node.local_UB = 0
        new_node.x_sol = copy.deepcopy(node.x_sol)
        new_node.x_int_sol = copy.deepcopy(node.x_int_sol)
        new_node.branch_var_list = []
        new_node.model = node.model.copy()
        new_node.cnt = node.cnt
        new_node.is_integer = node.is_integer

        return new_node

# Branch and Bound 
def Branch_and_bound(RLP):
    # initialise the initial node
    RLP.optimize()
    global_UB = RLP.ObjVal
    global_LB = -5000
    eps = 0.225
    incumbent_node = None
    Gap = None

    '''
        Branch and Bound starts
    '''
    # creat initial node
    Queue = []
    node = Node()
    node.local_LB = -2000
    node.local_UB = global_UB
    node.model = RLP.copy()
    node.model.setParam("OutputFlag", 0)
    node.cnt = 0
    Queue.append(node)

    cnt = 0
    Global_UB_change = []
    Global_LB_change = []

    while (len(Queue) > 0 and global_UB - global_LB > eps) :
        # select the current node
        current_node = Queue.pop()
        cnt += 1

        # solve the current model
        current_node.model.optimize()
        Solution_status = current_node.model.Status

        '''
        OPTIMAL = 2
        INFEASIBLE = 3
        UNBOUNDED = 5
        '''

        # check whether the current solution is integer and execute prune
        # step

        '''
            is_integer: mark whether the current solution is integer solution
            Is_Pruned: mark whether the current solution is pruned
        '''

        is_integer = True
        Is_Pruned = False
        
        if (Solution_status == 2):
            
            my_int_var_list=[]
            for i in range(2):
                for j in range(2):
                    my_int_var_list.append(current_node.model.getVarByName('rtsr_' + str(i)+ str(j)))
            for i in range(2):
                for j in range(2):
                    my_int_var_list.append(current_node.model.getVarByName('rtmf_' + str(i)+ str(j)))
            for i in range(2):
                for j in range(2):
                    my_int_var_list.append(current_node.model.getVarByName('rts_' + str(i)+ str(j)))
            for i in range(2):
                for j in range(2):
                    my_int_var_list.append(current_node.model.getVarByName('rta_' + str(i)+ str(j)))
            for i in range(2):
                my_int_var_list.append(current_node.model.getVarByName('qm_' + str(i)))
            for i in range(2):
                my_int_var_list.append(current_node.model.getVarByName('qs_' + str(i)))
            for i in range(2):
                my_int_var_list.append(current_node.model.getVarByName('ym_' + str(i)))
            for i in range(2):
                my_int_var_list.append(current_node.model.getVarByName('ys_' + str(i)))
            for i in range(2):
                my_int_var_list.append(current_node.model.getVarByName('dp_' + str(i)))
            for i in range(2):
                my_int_var_list.append(current_node.model.getVarByName('dm_' + str(i)))
            # 上面这一部分是筛选自己需要分支的变量,每一个'for' 循环都是用变量名将需要的变量从所有变量里提取出来
            
            for var in my_int_var_list:
                current_node.x_sol[var.varName] = var.x
                print (var.VarName, ' = ', var.x)
                
                current_node.x_int_sol[var.varName] = (int)(var.x) # round the solution to get an integer solution
                if (abs((int)(var.x) - var.x) >= eps):
                    is_integer = False
                    current_node.branch_var_list.append(var.VarName) # to record the candidate branch variables
            
            if (is_integer == True):
                # For integer solution rode, update the LB and UB
                current_node.is_integer = True
                current_node.local_LB = current_node.model.ObjVal
                current_node.local_UB = current_node.model.ObjVal
                if (current_node.local_LB > global_LB):
                    global_LB = current_node.local_LB
                    incumbent_node = Node.deepcopy_node(current_node)
            if (is_integer == False) :
            # For integer solution node, update the LB and UB also
                current_node.is_integer = False
                current_node.local_UB = current_node.model.ObjVal
                current_node.local_LB = -5000
                for var_name in current_node.x_int_sol.keys() :
                    var = current_node.model.getVarByName(var_name)
                    current_node.local_LB+=current_node.x_int_sol[var_name] * var.Obj

                if (current_node.local_LB > global_LB or (current_node.local_LB == global_LB and current_node.is_integer == True)):
                    global_LB = current_node.local_LB
                    incumbent_node = Node.deepcopy_node(current_node)
                    incumbent_node.local_LB = current_node.local_LB
                    incumbent_node.local_UB = current_node.local_UB
            '''
                PRUNE step
            '''
            # prune by optimility
            if (is_integer == True):
                Is_Pruned = True

            # prune by bound
            if (is_integer == False and current_node.local_UB < global_LB) :
                Is_Pruned = True

            Gap = round(100 * (global_UB - global_LB) / global_LB, 2)
            print('\n-------------- \n', cnt, '\t Gap = ', Gap,' %' )
        elif (Solution_status != 2):

            # the current node is infeasible or unbound
            is_integer = False
            '''
            PRUNE step
            '''
            # prune by infeasiblity
            Is_Pruned = True
            continue
        '''
            BRANCH step
        '''

        if (Is_Pruned == False) :
            # selecte the branch variable
            branch_var_name = current_node.branch_var_list[0]
            left_var_bound = (int)(current_node.x_sol[branch_var_name])
            right_var_bound = (int)(current_node.x_sol[branch_var_name]) + 1
            
            # creat two child nodes
            left_node = Node.deepcopy_node(current_node)
            right_node = Node.deepcopy_node(current_node)

            # creat left child node
            temp_var = left_node.model.getVarByName(branch_var_name)
            left_node.model.addConstr(temp_var <= left_var_bound, name='branch_left_' + str(cnt))
            left_node.model.setParam("OutputFlag", 0)
            left_node.model.update()
            cnt += 1
            left_node.cnt = cnt

            # creat right child node
            temp_var = right_node.model.getVarByName(branch_var_name)
            right_node.model.addConstr(temp_var >= right_var_bound, name = 'branch_right_' + str(cnt))
            right_node.model.setParam("OutputFlag",0)
            right_node.model.update()
            cnt += 1
            right_node.cnt=cnt
            
            Queue.append(left_node)
            Queue.append(right_node)

            # update the global UB, explor all the leaf nodes
            temp_global_UB = 0
            for node in Queue:
                node.model.optimize()
                if (node.model.status == 2):
                    if (node.model.ObjVal >= temp_global_UB):
                        temp_global_UB = node.model.ObjVal

            global_UB = temp_global_UB
            Global_UB_change.append(global_UB)
            Global_LB_change.append(global_LB)
            
      
    # all the nodes are explored, update the LB and UB
    global_UB = global_LB
    Gap = round (100 * (global_UB - global_LB) / global_LB, 2)
    Global_UB_change.append(global_UB)
    Global_LB_change.append(global_LB)

    print ('\n\n\n\n')
    print ('----------------------------------------')
    print ('      Branch and Bound terminates       ')
    print ('      Optimal solution found            ')
    print ('----------------------------------------')
    print ('\nFinal Gap= ', Gap, " %")
    print ('Optimal Solution:', incumbent_node.x_int_sol) 
    print ('Optimal ObJ:', global_LB) 

    return incumbent_node, Gap, Global_UB_change, Global_LB_change


# Solve the IP model by branch and bound
incumbent_node, Gap, Global_UB_change, Global_LB_change = Branch_and_bound(RLP)

# plot the results
# fig a plt.figure(r)
# plt.figure(figsize = (15, 10))
font_dict = {"family": 'Arial',   #"Kaiti"
            "style": "oblique",
            "weight": "normal",
            "color": "green",
            "size": 20
            }
plt.rcParams['figure.figsize'] = (12.0, 8.0)  #单位量inches
plt.rcParams["font.family"] = 'Arial' #"SimHei"
plt.rcParams["font.size"] = 16

x_cor = range(1, len(Global_LB_change) + 1)
plt.plot(x_cor, Global_LB_change, label = 'LB')
plt.plot(x_cor, Global_UB_change, label = 'UB')
plt.legend()
plt.xlabel('Iteration', fontdict=font_dict)
plt.ylabel('Bounds update', fontdict=font_dict)
plt.title('Bounds update during branch and bound procedure \n', fontsize = 23)
plt.savefig('Bound_updates.eps')
plt.show()
