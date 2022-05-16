
from distutils.util import execute
from pkgutil import extend_path
from cv2 import ml_ANN_MLP, solve
from pandas import option_context
from pyomo.environ import *
import numpy as np


def lowbit(x):
    return (-x) & x


N = 10
L = np.random.randint(low=N, size=N)
R = np.random.randint(low=N, size=N)
# path = "/Users/nacn/Graduation-Design/ipopt-osx/ipopt"
# Q = [[0 for i in range(N)] for j in range(N)]
B = [[0 for i in range(N)] for j in range(N)]
M = [[0 for i in range(N)] for j in range(N)]

for i in range(1, N+1):
    if(L[i-1] > R[i-1]):
        L[i-1], R[i-1] = R[i-1], L[i-1]
    j = R[i-1]+1

    while j > 0:
        B[i-1][j-1] += 1.0
        j -= lowbit(j)

    j = L[i-1]
    while j > 0:
        B[i-1][j-1] -= 1.0
        j -= lowbit(j)

for i in range(1, N+1):
    j = i
    while j <= N:
        M[j-1][i-1] = 1.0
        j += lowbit(j)
B = np.array(B)
M = np.array(M)
print(B)
print(M)
print(B.dot(M))

model = ConcreteModel()

model.I = Set(initialize=[i for i in range(N)])
model.J = Set(initialize=[i for i in range(N)])

model.lamda = Var(model.I, within=PositiveReals)
# model.l = Var(model.I, model.J)


def obj_rule(model):
    sum = 0
    for i in range(N):
        for j in range(N):
            sum += 1/(model.lamda[j]*model.lamda[j])*B[i][j]*B[i][j]
    return sum


def c2_rule(model,  j):
    sum = 0
    for i in range(N):
        sum += model.lamda[i]*M[i][j]
    return sum <= 1


# model.c1 = Constraint(model.I, model.J, rule=c1_rule)
model.c2 = Constraint(model.J, rule=c2_rule)
model.obj = Objective(rule=obj_rule, sense=minimize)


opt = SolverFactory('ipopt')
# opt.options['max_iter'] = 50
model.pprint()
opt.solve(model, tee=True)
for i in range(N):
    for j in range(N):
        B[i][j] *= 1/value(model.lamda[j])
for i in range(N):
    for j in range(N):
        M[i][j] *= value(model.lamda[i])
# print(B)
# print(M)
# print(B.dot(M))
ansL = np.array([value(model.lamda[i])for i in model.I])
# # ansB = ansB.reshape((N, N))
print(ansL)
# # for j in range(N):
# #     sum = 0
# #     for i in range(N):
# #         sum += abs(value(model.l[i, j]))
# #     print(sum)
# # ansL = np.array([value(model.l[i, j])for i in model.I
# #                  for j in model.J])
# # ansL = ansL.reshape((N, N))
# # print(ansB)
# # print(ansL)
# # print(np.matmul(ansB, ansL))

print(model.obj())
