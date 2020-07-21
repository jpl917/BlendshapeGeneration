import numpy as np
from config import Config
basepath = 'E:\Tester_1'
objFilePath = 'E:\Tester_1\Blendshape\shape_0.obj'
posepath = basepath + '\TrainingPose'
bspath = basepath + '\Blendshape'

face_vertex = [] #三角面片与顶点的对应关系
edge_vpair = [] #几何上可以理解为边与顶点的对应关系，在求解localframe时边即为基本的未知参数
face_edge_index = [] #三角面片与lf中边的对应关系，一对二
pose_vlist = [] #各pose包含的顶点坐标，pose_vlist[0][1]对应第一个pose的第二个顶点
pose_lflist = [] #各pose的localframe信息，pose_lflist[0][1]对应第一个pose的第二个三角面片的localframe(3*3)


def localframe(tri,pi): #求解lf的函数
    data = np.array(pose_vlist[pi])
    a = data[tri[2]][0] - data[tri[0]][0]
    b = data[tri[1]][0] - data[tri[0]][0]
    n = np.cross(a,b)
    lf = np.array([a,b,n])
    return lf


def g_0i(mb0, mbi):
    return g_ab(mb0, mbi) + g_ab(mb0, mb0)
def g_ab(ma, mb):
    A = np.mat(ma)
    B = np.mat(mb)
    A = A.astype(np.float).I
    B = B.astype(np.float)
    return np.dot(B,A)

i = 0
for pi in range(0,Config.POSENUM):
    posename = '\pose_' + str(pi) + '.obj'
    with open(posepath + posename) as file:
        points = []
        vt = []
        face_vertex = []
        fvt = []
        while 1:
            line = file.readline()
            sfv1 = []
            sfv2 = []
            sfvt1 = []
            sfvt2 = []
            if not line:
                break
            strs = line.split(" ")
            if strs[0] == "v":
                points.append([(float(strs[1]), float(strs[2]), float(strs[3]))])
            if strs[0] == "vt":
                vt.append([float(strs[1]), float(strs[2])])
            if strs[0] == "f":
                for i in range(1, 5):  # 将四边形切分为三角面片
                    dum = strs[i].split("/")
                    if (i != 4):
                        sfv1.append(int(dum[0]))
                        sfvt1.append(int(dum[1]))
                    if (i != 2):
                        sfv2.append(int(dum[0]))
                        sfvt2.append(int(dum[1]))
                if(pi == 0):
                    face_vertex.append(sfv1)
                    face_vertex.append(sfv2)
    '''''''''
    if (pi == 0):  # 各model的拓扑逻辑是一致的，因此如面点对应关系、localframe序列等信息只需要统计一次
        print("统计")
        face_edge_index = []
        for face in face_vertex:
            face_edge_index.append([0, 0])
        i = 0
        fi = 0
        for face in face_vertex:
            edge1 = [face[2], face[0]]
            edge2 = [face[1], face[0]]  # lf的两个对应边
            if (edge1 not in edge_vpair):
                edge_vpair.append(edge1)
                face_edge_index[fi][0] = i
                i = i + 1
            else:
                face_edge_index[fi][0] = edge_vpair.index(edge1)
            if (edge2 not in edge_vpair):
                edge_vpair.append(edge2)
                face_edge_index[fi][1] = i
                i = i + 1
                #
            else:
                face_edge_index[fi][1] = edge_vpair.index(edge2)
            fi = fi + 1
            print(fi)
        edge_vpair = np.array(edge_vpair)
        face_edge_index = np.array(face_edge_index)
        np.save('ep.npy', edge_vpair)
        np.save('fe.npy', face_edge_index)
        '''''''''
    edge_vpair = np.load('ep.npy')
    face_edge_index = np.load('fe.npy')
    pose_vlist.append(points)
    localframes = []
    for f in face_vertex:
        res = localframe(f, pi)
        localframes.append(res)
        print(res)
    pose_lflist.append(localframes)
    print(len(pose_lflist))
points = np.array(points)


