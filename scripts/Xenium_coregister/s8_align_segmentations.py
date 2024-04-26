#!/usr/bin/env python3

# THIS FILE WILL ALIGN A SEGMENTATION WITH A TARGET PNG IMAGE. 
# LANDMARKS ARE REQUIRED: ASK JUNXIANG ABOUT THIS

# USAGE FOR JUNXIANG ALIGNED XENIUM DATA (DEFAULT):
# ./s8_align_segmentation.py -s segmentation_file_to_align.tif
# OR FOR ALTERNATE TARGET ./s8_align_segmentation.py -t target_image.png -s segmentation_file.tif 



import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.append("./STalign/")
import torch
# import STalign from STalign directory
import STalign
plt.rcParams["figure.figsize"] = (12,10)
from STalign import STalign


# Default Target is H&E staining image
image_file = "/projectnb/rd-spat/HOME/junxiang/Image_Registration/Xenium_register/Images/Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image.png"
fname = '/projectnb/rd-spat/HOME/junxiang/Image_Registration/Xenium_register/Images/IF_DAPI_point_rep1.csv'
seg_fname = '/projectnb/rd-spat/HOME/junxiang/Image_Registration/Xenium_register/Segmentation/star_dist_seg.csv'
if_dapi_landmarks = '/projectnb/rd-spat/HOME/junxiang/Image_Registration/Aligned_Xe_rep1/points/IF_DAPI_point_points.npy'
he_landmarks = '/projectnb/rd-spat/HOME/junxiang/Image_Registration/Aligned_Xe_rep1/points/Xenium_HE_points.npy'

parser = argparse.ArgumentParser("Align IF-based Segmentation Data using STalign.")
parser.add_argument("-t",
                    "--target",
                    type=str,
                    default = image_file)
parser.add_argument("-if",
                    "--if_image",
                    type=str, 
                    default = fname)
parser.add_argument("-s",
                    "--segmentation", 
                    type=str, 
                    default = seg_fname)
parser.add_argument("-o",
                    "--output",
                    type=str,
                    default = "aligned_segmentation.csv")
parser.add_argument("-i",
                    "--i_landmark_points",
                    type=str,
                    default = if_dapi_landmarks)
parser.add_argument("-j",
                    "--j_landmark_points",
                    type=str,
                    default = he_landmarks)
parser.add_argument("-p",
                    "--save_and_plot_image",
                    type=str,
                    default = "True")


def read_landmark_files(if_dapi_landmarks: str, he_landmarks: str):
    ## Read landmarks
    
    # read from file
    pointsIlist = np.load(if_dapi_landmarks, allow_pickle=True).tolist()
    pointsJlist = np.load(he_landmarks, allow_pickle=True).tolist()
    
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
    
    return pointsI, pointsJ



def affine_transform_segmentation(df: pd.DataFrame,
                                  L: np.ndarray,
                                  T: np.ndarray,
                                  V: np.ndarray ):
    
    xSeg = np.array(df['x'])
    ySeg = np.array(df['y'])
    
    # note points are as y,x
    affine = np.dot(np.linalg.inv(L), [ySeg - T[0], xSeg - T[1]])
    x_seg_affine = V.shape[0] - affine[0,:]
    y_seg_affine = affine[1,:]
    
    return x_seg_affine, y_seg_affine


if __name__ == "__main__":
    args = parser.parse_args()
    
    image_file = args.target
    fname = args.if_image
    seg_fname = args.segmentation
    output_file_name = args.output
    if_dapi_landmarks = args.i_landmark_points
    he_landmarks = args.j_landmark_points
    save_and_plot = args.save_and_plot_image
    
    if (save_and_plot == "True"): save_and_plot = True
    else: save_and_plot = False
    
    print(f"Beginning segmentation alignment for: \n{seg_fname}.")
    # READ IN IMAGE
    V = plt.imread(image_file)

    # Normalize image
    Jnorm = STalign.normalize(V)
    
    #Transpose Image to a np array
    J = Jnorm.transpose(2,0,1)
    #print(J.shape)
    
    YJ = np.array(range(J.shape[1]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
    XJ = np.array(range(J.shape[2]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
    extentJ = STalign.extent_from_x((YJ,XJ))
    
    # IF Data to be aligned
    
    df_if = pd.read_csv(fname)
    #filter by intensity
    df_if = df_if[df_if['x'] > 0.1]
    
    # get cell centroid coordinates
    xD = np.array(df_if['i'])
    yD = np.array(df_if['j'])
    
    # Rasterize
    XD,YD,D,fig = STalign.rasterize(xD, yD, dx=1)
    
    # Process Raster image to be a 3*N*M npy array,then normalize
    D = np.vstack((D, D, D)) # make into 3xNxM
    # normalize
    D = STalign.normalize(D)
    
    pointsI, pointsJ = read_landmark_files(if_dapi_landmarks, he_landmarks)
    print("Preparing Torch Default device.")
    # set device for building tensors
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
    
    #Point transfer
    L,T = STalign.L_T_from_points(pointsJ,pointsI)
    
    # Read in segmentation data
    df = pd.read_csv(seg_fname)
    
    x_seg_affine, y_seg_affine = affine_transform_segmentation(df, L, T, V)
    
    if(save_and_plot):
        fig,ax = plt.subplots()
        ax.scatter(y_seg_affine,x_seg_affine,s=0.01,alpha=0.2)
        plt.savefig('aligned_segmentation.png',dpi = 800)
    
    df_seg = pd.DataFrame({"aligned_x": y_seg_affine,"aligned_y": x_seg_affine})
    seg_poly = pd.concat([df, df_seg], axis=1)
    seg_poly.to_csv(output_file_name)
    
    print(f"Saved aligned segmentation data as {output_file_name}.")
