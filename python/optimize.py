import scipy
import scipy.optimize as op
import numpy as np
from cvxopt import solvers
from config import Config

#求解lf的优化函数
def solve_localframe(restpose, pose_list, model_list, weight, edge_vpair, face_edge_index):
    lf_num = len(edge_vpair)
    n = lf_num * 3 * Config.BLENDSHAPENUM
    m = lf_num * 3 * Config.BLENDSHAPENUM
    A = np.zeros(m, n)  # Ax=b
    B = np.zeros(n, 1)
    for fi in range(0, len(restpose)):
        for i1 in range(0, Config.BLENDSHAPENUM):
            for t in range(0, 2):
                w = (1 + np.linalg.norm(model_list[i1][fi])) / (Config.Ka + np.linalg.norm(model_list[i1][fi]))
                w = w ** Config.Theta

                A[(lf_num * i1 + face_edge_index[fi][t]) * 3][(lf_num * i1 + face_edge_index[fi][t]) * 3] += 2 * w
                B[(lf_num * i1 + face_edge_index[fi][t]) * 3] += 2 * model_list[i1][fi][t][0] * w

                A[(lf_num * i1 + face_edge_index[fi][t]) * 3 + 1][(lf_num * i1 + face_edge_index[fi][t]) * 3 + 1] += 2 * w
                B[(lf_num * i1 + face_edge_index[fi][t]) * 3 + 1] += 2 * model_list[i1][fi][t][1] * w

                A[(lf_num * i1 + face_edge_index[fi][t]) * 3 + 2][(lf_num * i1 + face_edge_index[fi][t]) * 3 + 2] += 2 * w
                B[(lf_num * i1 + face_edge_index[fi][t]) * 3 + 2] += 2 * model_list[i1][fi][t][2] * w
            for j in range(0, len(pose_list)):
                for t in range(0, 2):
                    B[(lf_num * i1 + face_edge_index[fi][t]) * 3] += (restpose[fi][t][0] - pose_list[j][fi][t][0]) * weight[j][i1] * 2
                    B[(lf_num * i1 + face_edge_index[fi][t]) * 3 + 1] += (restpose[fi][t][1] - pose_list[j][fi][t][1]) * weight[j][i1] * 2
                    B[(lf_num * i1 + face_edge_index[fi][t]) * 3 + 2] += (restpose[fi][t][1] - pose_list[j][fi][t][1]) * weight[j][i1] * 2
                for i2 in range(i1,Config.BLENDSHAPENUM):
                    for t in range(0, 2):
                        A[(lf_num * i1 + face_edge_index[fi][t]) * 3][(lf_num * i2 + face_edge_index[fi][t]) * 3] += 2 * weight[j][i2] * weight[j][i1]
                        A[(lf_num * i1 + face_edge_index[fi][t]) * 3 + 1][(lf_num * i2 + face_edge_index[fi][t]) * 3 + 1] += 2 * weight[j][i2] * weight[j][i1]
                        A[(lf_num * i1 + face_edge_index[fi][t]) * 3 + 2][(lf_num * i2 + face_edge_index[fi][t]) * 3 + 2] += 2 * weight[j][i2] * weight[j][i1]

    result =  scipy.linalg.solve(A, B)
    solved_localframe = np.zeros(lf_num,3)
    output = []
    count = 0
    index = 0
    unit = np.zeros(1,3)
    for item in result:
        unit[count] = item
        count += 1
        if count == 3:
            solved_localframe[index] = unit
            index += 1
        if index == lf_num:
            output.append(solved_localframe)
            solved_localframe = np.zeros(lf_num,3)
            index = 0
    #bsnum * lfnum * 3
    def E_fit():
        sum = 0
        target = [0, 0, 0]
        c = [0, 0, 0]
        for j in range(0, len(pose_list)):
            for fi in range(0, len(restpose)):
                for t in range(0, 2):
                    c = restpose[fi][t] - pose_list[j][fi][t]
                    for i in range(0, Config.BLENDSHAPENUM):
                        target += output[i][face_edge_index[fi][t]] * weight[j][i]
        res = c - target
        sum += (res[0] ** 2 + res[1] ** 2 + res[2] ** 2)
        return sum

    def E_reg():
        sum = 0
        for fi in range(0, len(restpose)):
            for i in range(0, Config.BLENDSHAPENUM):
                w = (1 + np.linalg.norm(model_list[i][fi])) / (Config.Ka + np.linalg.norm(model_list[i][fi]))
                w = w ** Config.Theta
                for t in range(0, 2):
                    target = output[i][face_edge_index[fi][t]] - model_list[i][fi][t]
                    res = target[0] ** 2 + target[1] ** 2 + target[2] ** 2
                    sum = sum + w * res
        return sum
    print("Energy: " + str(E_fit()+E_reg()))

#求解权重的函数
def solve_weight(pose_list,bs_model,blendshape_list,pre_weight):
    weightnum = Config.BLENDSHAPENUM * len(pose_list)
    bsnum = Config.BLENDSHAPENUM
    P = np.zeros(weightnum,weightnum)
    Q = np.zeros(weightnum,1)
    G = np.zeros(weightnum,weightnum) - np.eye(weightnum,weightnum)
    H = np.zeros(weightnum,1)
    A = np.ones(weightnum,1)
    B = [1.0]
    for j in range(0,len(pose_list)):
        for k in range(0,len(pose_list[0])):
            c = pose_list[j][k] - bs_model[k]
            for i1 in range(0,Config.BLENDSHAPENUM):
                Q[j * bsnum + i1] += 2 * c[0] * blendshape_list[i1][k][0]
                Q[j * bsnum + i1] += 2 * c[1] * blendshape_list[i1][k][1]
                Q[j * bsnum + i1] += 2 * c[2] * blendshape_list[i1][k][2]
                for i2 in range(i1,Config.BLENDSHAPENUM):
                    P[j * bsnum + i1][j * bsnum + i2] += blendshape_list[i1][k][0] * blendshape_list[i2][k][0] * 2
                    P[j * bsnum + i1][j * bsnum + i2] += blendshape_list[i1][k][1] * blendshape_list[i2][k][1] * 2
                    P[j * bsnum + i1][j * bsnum + i2] += blendshape_list[i1][k][2] * blendshape_list[i2][k][2] * 2
                    if(i1 != i2):
                        P[j * bsnum + i2][j * bsnum + i1] += blendshape_list[i1][k][0] * blendshape_list[i2][k][0] * 2
                        P[j * bsnum + i2][j * bsnum + i1] += blendshape_list[i1][k][1] * blendshape_list[i2][k][1] * 2
                        P[j * bsnum + i2][j * bsnum + i1] += blendshape_list[i1][k][2] * blendshape_list[i2][k][2] * 2

    for j in range(0,len(pose_list)):
        for i in range(0, Config.BLENDSHAPENUM):
            P[j * bsnum + i][j * bsnum + i] += 2 * Config.Nju
            Q[j * bsnum + i] += -2 * Config.Nju * pre_weight[j][i]

    result = solvers.qp(P, Q, G, H, A, B)

    def E_weight():
        N = pose_list[0].length()
        sum1 = 0.0
        for k in range(N):
            for j in range(len(pose_list)):
                target = [0, 0, 0]
                for i in range(len(blendshape_list)):
                    target += result[j * bsnum + i] * blendshape_list[j][i][k]
                target = pose_list[j][k] - bs_model[k] - target
                sum1 += (target[0] ** 2 + target[1] ** 2 + target[2] ** 2)
        sum2 = 0.0
        for j in range(len(pose_list)):
            for i in range(weightnum):
                t = result[j * bsnum + i] - pre_weight[j][i]
                sum2 += t ** 2
        sum2 *= Config.Nju
        return sum1 + sum2

    print("Weight energy: " + str(E_weight()))
    return result

#从localframe重构顶点坐标的函数
"""""""""""
def recover_mesh_by_localframe(lf,lf_v,x):
    
    def func(x,lf,lf_v):
        Error = []
        i = 0
        ii = 0
        for face in lf_v:
            v1 = face[0] * 3
            v2 = face[1] * 3
            cost = ((x[v1] - x[v2]) - lf[i][ii][0]) ** 2 + ((x[v1+1] - x[v2+1]) - lf[i][ii][1]) ** 2 + ((x[v1+2] - x[v2+2]) - lf[i][ii][2]) ** 2
            Error.append(cost)
            if ii == 0:
                ii = 1
            if ii == 1:
                i = i+1
                ii = 0
        return Error
    res = op.leastsq(func,x,args=(lf,lf_v))
    i = 0
    result = []
    unit = []
    for item in res:
        unit.append(item)
        if i == 2:
            result.append(unit)
            unit.clear()
        else:
            i = i + 1
    
    return result
"""""""""""