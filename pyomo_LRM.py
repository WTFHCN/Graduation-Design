
from pyomo.environ import *
import numpy as np



def calc(Q, N, M):

    Q = [[0 for i in range(N)] for j in range(M)]
    for i in range(N):
        for j in range(i+1):
            Q[i][j] = 1

    model = ConcreteModel()

    model.I = Set(initialize=[i for i in range(N)])
    model.J = Set(initialize=[i for i in range(M)])

    model.b = Var(model.I, model.J)
    model.l = Var(model.I, model.J)

    def obj_rule(model):
        sum = 0
        for i in range(N):
            for j in range(M):
                sum += model.b[i, j]*model.b[i, j]
        return sum

    def c1_rule(model, i, j):
        sum = 0
        for k in range(M):
            sum += model.b[i, k]*model.l[k, j]
        return sum == Q[i][j]

    def c2_rule(model,  j):
        sum = 0
        for i in range(N):
            sum += abs(model.l[i, j])
        return sum <= 1

    model.c1 = Constraint(model.I, model.J, rule=c1_rule)
    model.c2 = Constraint(model.J, rule=c2_rule)
    model.obj = Objective(rule=obj_rule, sense=minimize)

    opt = SolverFactory('ipopt')
    opt.options['max_iter'] = 50
    model.pprint()
    opt.solve(model, tee=True)
    ansB = np.array([value(model.b[i, j])for i in model.I
                     for j in model.J])
    ansB = ansB.reshape((N, N))

    for j in range(N):
        sum = 0
        for i in range(M):
            sum += abs(value(model.l[i, j]))
        print(sum)
    ansL = np.array([value(model.l[i, j])for i in model.I
                     for j in model.J])
    ansL = ansL.reshape((N, N))
    print(ansB)
    print(ansL)
    # print(np.matmul(ansB, ansL))

    print(model.obj())
