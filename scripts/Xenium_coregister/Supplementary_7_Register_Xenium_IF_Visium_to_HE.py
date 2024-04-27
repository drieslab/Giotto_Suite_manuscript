
# import dependencies
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.append("./STalign/")
import torch
# import STalign from STalign directory
import STalign

import argparse

parser.add_argument('input_file', type=str, help='The input file name')
parser.add_argument('output_file', type=str, help='The output file name')
    
args = parser.parse_args()

# make plots bigger
plt.rcParams["figure.figsize"] = (12,10)

from STalign import STalign

# Target is H&E staining image
image_file = "./Images/Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image_compressed.png"
V = plt.imread(image_file)

# Normalize image
Jnorm = STalign.normalize(V)

fig,ax = plt.subplots()

#Transpose Image to a np array
J = Jnorm.transpose(2,0,1)

YJ = np.array(range(J.shape[1]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
XJ = np.array(range(J.shape[2]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
extentJ = STalign.extent_from_x((YJ,XJ))
# IF DAPI data to be aligned
# Note that IF data needs to be converted from image to points first
fname = './Xenium_register/Images/IF_DAPI_point_rep1.csv'
df = pd.read_csv(fname)
#filter by intensity
df = df[df['x'] > 0.08]

# get cell centroid coordinates
xD = np.array(df['i'])
yD = np.array(df['j'])
# rasterize and plot
XD,YD,D,fig = STalign.rasterize(xD, yD, dx=5)
D = np.vstack((D, D, D)) # make into 3xNxM
D = STalign.normalize(D)
# get extent of images
extentJ = STalign.extent_from_x((YJ,XJ))
extentD = STalign.extent_from_x((YD,XD))
## Read landmarks and plot

# read from file
pointsIlist = np.load('./IF_DAPI_point_points.npy', allow_pickle=True).tolist()
pointsJlist = np.load('./Xenium_HE_points.npy', allow_pickle=True).tolist()

# convert to array
pointsI = []
pointsJ = []

for i in pointsIlist.keys():
    for j in range(len(pointsIlist[i])):
        pointsI.append([pointsIlist[i][j][1], pointsIlist[i][j][0]])
for i in pointsJlist.keys():
    for j in range(len(pointsJlist[i])):
        pointsJ.append([pointsJlist[i][j][1], pointsJlist[i][j][0]])

pointsI = np.array(pointsI) 
pointsJ = np.array(pointsJ) 

#For original image, times the scale factor
scale_factor = 1
pointsI = pointsI * scale_factor
# Using a dictionary comprehension with tuple comprehension to multiply each element by 10
pointsIlist = {k: [(x*scale_factor, y*scale_factor) for x, y in v] for k, v in pointsIlist.items()}
if torch.cuda.is_available():
    torch.set_default_device('cuda:0')
else:
    torch.set_default_device('cpu')
# compute initial affine transformation from points
L,T = STalign.L_T_from_points(pointsI,pointsJ)
A = STalign.to_A(torch.tensor(L),torch.tensor(T))
# We can show the results of the simple affine transformation.

# compute initial affine transformation from points
AD= STalign.transform_image_atlas_with_A(A, [YD,XD], D, [YJ,XJ])
#The results seems no velocity needed

#Point transfer
L,T = STalign.L_T_from_points(pointsJ,pointsI)
# note points are as y,x
affine = np.dot(np.linalg.inv(L), [yD- T[0], xD- T[1]])
xDaffine = affine[0,:]
yDaffine = affine[1,:]

# Apply to all IF image based segmentation
# Star Dist segmentation to be aligned
fname = args.input_file

df = pd.read_csv(fname)
xStarDist = np.array(df['x'])
yStarDist = np.array(df['y'])

affine = np.dot(np.linalg.inv(L), [yStarDist - T[0], xStarDist - T[1]])
x_StarDist_affine = V.shape[0] - affine[0,:]
y_StarDist_affine = affine[1,:]

df_StarDist = pd.DataFrame({"aligned_x": y_StarDist_affine,"aligned_y": x_StarDist_affine})
StarDist_poly = pd.concat([df, df_StarDist], axis=1)
StarDist_poly.to_csv(args.output_file,compression='gzip')
