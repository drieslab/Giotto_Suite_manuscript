## import dependencies
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import pandas as pd

import sys
sys.path.append("./STalign/")

import torch
import plotly
import requests

# make plots bigger
plt.rcParams["figure.figsize"] = (12,10)


from STalign import STalign

# Target is H&E staining image
image_file = "../Visium/spatial/tissue_hires_image.png"
V = plt.imread(image_file)

# Normalize image
Jnorm = STalign.normalize(V)

#Transpose Image to a np array
J = Jnorm.transpose(2,0,1)
YJ = np.array(range(J.shape[1]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
XJ = np.array(range(J.shape[2]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
extentJ = STalign.extent_from_x((YJ,XJ))

# Single cell data to be aligned
fname = './Xenium_register/outs_rep1/cells.csv.gz'
df = pd.read_csv(fname)
# get cell centroid coordinates
# flip x and y to make it align
xI = np.array(df['y_centroid'])
yI = np.array(df['x_centroid'])

# rasterize and plot
XI,YI,I,fig = STalign.rasterize(xI, yI, dx=5)
I = np.vstack((I, I, I)) # make into 3xNxM
# normalize
# Normalize image
I = STalign.normalize(I)
extentI_tx = STalign.extent_from_x((YI,XI))

## Read landmarks and plot

# read from file
pointsIlist_tx = np.load('./xenium_rep1_to_visium/points/Xenium_to_VisiumHE_points.npy', allow_pickle=True).tolist()
pointsJlist_tx = np.load('./xenium_rep1_to_visium/points/Visium_points.npy', allow_pickle=True).tolist()

# convert to array
pointsI_tx = []
pointsJ_tx = []

for i in pointsIlist_tx.keys():
    for j in range(len(pointsIlist_tx[i])):
        pointsI_tx.append([pointsIlist_tx[i][j][1], pointsIlist_tx[i][j][0]])
for i in pointsJlist_tx.keys():
    for j in range(len(pointsJlist_tx[i])):
        pointsJ_tx.append([pointsJlist_tx[i][j][1], pointsJlist_tx[i][j][0]])

pointsI_tx = np.array(pointsI_tx)
pointsJ_tx = np.array(pointsJ_tx)
# compute initial affine transformation and show points
L_tx,T_tx = STalign.L_T_from_points(pointsJ_tx,pointsI_tx)
# note points are as y,x
affine = np.dot(np.linalg.inv(L_tx), [yI- T_tx[0] ,xI - T_tx[1]])
xIaffine = V.shape[0] - affine[0,:]
yIaffine = affine[1,:]

L_poly,T_poly = STalign.L_T_from_points(pointsI_tx,pointsJ_tx)
affine = np.dot(np.linalg.inv(L_poly), [Y_poly - T_poly[0], X_poly - T_poly[1]])
# Keep compute initial affine transformation for LDDMM(image fields)
# **Note the change of landmark order!!!
L_tx,T_tx = STalign.L_T_from_points(pointsI_tx,pointsJ_tx)
A_tx = STalign.to_A(torch.tensor(L_tx),torch.tensor(T_tx))
# compute initial affine transformation from points
AI_tx= STalign.transform_image_atlas_with_A(A_tx, [YI,XI], I, [YJ,XJ])

%%time
# set device for building tensors
if torch.cuda.is_available():
    torch.set_default_device('cuda:0')
else:
    torch.set_default_device('cpu')

# run LDDMM
# specify device (default device for STalign.LDDMM is cpu)
if torch.cuda.is_available():
    device = 'cuda:0'
else:
    device = 'cpu'

# keep all other parameters default
params = {'L':L_tx,'T':T_tx,
          'niter': 200,
          'pointsI': pointsI_tx,
          'pointsJ': pointsJ_tx,
          'device': device,
          'sigmaP': 2e-1,
          'sigmaM': 0.18,
          'sigmaB': 0.18,
          'sigmaA': 0.18,
          'diffeo_start' : 100,
          'epL': 5e-11,
          'epT': 5e-4,
          'epV': 5e1
          }

A_tx,v_tx,xv_tx = STalign.LDDMM([YI,XI],I,[YJ,XJ],J,**params)

# apply transform
phii_tx = STalign.build_transform(xv_tx,v_tx,A_tx,XJ=[YJ,XJ],direction='b')
phiI_tx = STalign.transform_image_atlas_to_target(xv_tx,v_tx,A_tx,[YI,XI],I,[YJ,XJ])
phiipointsI_tx = STalign.transform_points_atlas_to_target(xv_tx,v_tx,A_tx,pointsI_tx)

# apply transform to original points
tpointsI_tx= STalign.transform_points_atlas_to_target(xv_tx,v_tx,A_tx, np.stack([yI, xI], 1))

# switch from row column coordinates (y,x) to (x,y)
xI_LDDMM = tpointsI_tx[:,1]
yI_LDDMM = V.shape[0] - tpointsI_tx[:,0]

if phiipointsI_tx.is_cuda:
    landmarks_tx_x = phiipointsI_tx[:,1].cpu().detach()
    landmarks_tx_y = V.shape[0] - phiipointsI_tx[:,0].cpu().detach()
else:
    landmarks_tx_x = phiipointsI_tx[:,1].cpu().detach()
    landmarks_tx_y = V.shape[0] - phiipointsI_tx[:,0].cpu().detach()

landmarks_HE_x = pointsJ_tx[:,1]
landmarks_HE_y = V.shape[0] - pointsJ_tx[:,0]
df3 = pd.DataFrame(

        {

            "aligned_x": xI_LDDMM.cpu(),

            "aligned_y": yI_LDDMM.cpu(),

        }
)
results = pd.concat([df, df3], axis=1)
results.to_csv('./Xenium_Rep1_STalign_to_VisiumHE.csv')


