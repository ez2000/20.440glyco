#########################
### PACKAGE SCRIPT ###
#########################

#-----------------------------------------------------------------
# 1) Install packages
#-----------------------------------------------------------------

### Install and load packages quietly (console doesn't get messy) 
library(remotes)
library(fields)

##from Bojar lab github
library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(DropletUtils)
library(Matrix)
library(glmGamPoi)
library(ProjecTILs) #download package from https://github.com/carmonalab/ProjecTILs
library(Cairo)
library(ggplot2)
library(RColorBrewer)
library(scater)
library(rgdal)

suppressMessages(require(Seurat))
suppressMessages(require(dplyr))
suppressMessages(require(patchwork))
suppressMessages(require(Matrix))
