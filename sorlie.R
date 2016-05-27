# Copyright 2016 The Hyve BV.
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version, along with the following terms:
#
#   1. You may convey a work based on this program in accordance with
#      section 5, provided that you retain the above notices.
#   2. You may convey verbatim copies of this program code as you receive
#      it, in any medium, provided that you retain the above notices.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <http://www.gnu.org/licenses/>.

# This file contains an example of using TranSMART RClient to fetch study data.

#############################################################################
#                             Introduction                                  #
#############################################################################
# Welcome to the Rclient demo for the Open Data workshop. For the workshop
# we will use the Sorlie(2003) study to demonstrate how to retrieve the data
# from a tranSMART instance and generate a heatmap basd on the actual paper.
# 
# This server with Rstudio, Rscripts and Rclient will stay available until
# Friday 3rd of June. If you want to redo the analysis or have a closer look
# at the code please feel free do download and use. Note: You will have to
# install the tranSMART RClient on your own machine.
# For transmartRClient installation see https://github.com/thehyve/RInterface
#
# For any questions please contact us at office@thehyve.nl or support@thehyve.nl

#############################################################################
#                               Workshop                                    #
#############################################################################
# Example steps to authenticate with, connect to, and retrieve data from tranSMART
# load package
require("transmartRClient")

# Connect to the Hyve's demo transmart instance with the connectToTransmart function.
#   You can log in as a guest user with username "user" and password "user". 
connectToTransmart("http://transmart-demo.thehyve.net/transmart")

# retrieve a list of the available studies in the database:
studies <- getStudies()
print(studies)
View(studies)

# to access the studies programmatically use for example: 
# study<-studies$id[1]
# for the examples below we will use GSE4382, which is called "Sorlie(2003) GSE4382" in transmart.
# This also shows that a study name is not the same as a study ID, all the data retrieval is done
# study ID rather than the name. You can use View() to have a closer look at the studies in transmart.
study <- "GSE4382_ETRIKS"  

# Retrieve Clinical Data
allObservations <- getObservations(study, as.data.frame = T)

#the clinical data is stored in three tables: 
# - an observation table with the clinical measurements
# - a subjectInfo table which contains information about the subjects/patients, such as their age or sex
# - and a table with information about the variables/concepts and where they can be 
#    found in the data tree  of tranSMART (see the panel "Navigate terms" in the "Analyze" tab of the 
#    tranSMART web app: http://transmart-demo.thehyve.net/transmart/datasetExplorer/index ) 
summary(allObservations)
allObservations$observations[1:12,1:7] #hint: if you are using RStudio you can also use the function "View" to see the data in a more user-friendly table, e.g.: View(allObservations$observations)
allObservations$subjectInfo[1:12,]
allObservations$conceptInfo[1:10,]

# retrieve information about the concepts that are present for this study
concepts <- getConcepts(study)
View(concepts[c("name","type","fullName")])

# You can use the concept names or the concept links to only retrieve data for a subset of the variables:

# retrieve observations for the first concept which concept name contains "Tumor Subtypes"
observations <- getObservations(study, concept.match = "Tumor Subtypes", as.data.frame = T)
View(observations$observations)

# retrieve observations belonging to the first two concepts by using the api.links contained in the getConcepts-result
observations <- getObservations(study, concept.links = concepts$api.link.self.href[c(1,2)])
View(observations$observations[1:10,])


# if a concept contains high dimensional data, use the following command to obtain this data. 
# NB: you will be told that you need to select one of the listed projections. For the full command, 
# see next line below. 
getHighdimData(study.name = study, concept.match = "Breast")

# As noted above, you will be told that one of the listed projections needs to be selected. 
# The following will return the actual data, log transformed.
dataDownloaded <- getHighdimData(study.name = study, concept.match = "Breast", projection = "log_intensity")

# getHighDimData returns a list containing two objects: 
# a data.frame containing the data, and a hash which maps probe names (labels) to Biomarker (e.g. gene) names
summary(dataDownloaded)

#View the data
View(dataDownloaded$data[1:10,1:10])

#select gene expression data
expression_data<-dataDownloaded$data[,-c(1:6)]
dim(expression_data)
rownames(expression_data)<-dataDownloaded$data$patientId
View(expression_data[1:10,1:5])

#Make a heatmap
# Please note that if the dimensions of the expression_data table are large,
# it may take some time to generate the heatmap. You may want to create a subset
# of the data first, which we will do now for the Sorlie(2003) example.

# Create a dataframe where it is clear which patients are part of what group,
# remove the patients that do not have a subgroup ("None")
# For this step you use the downloaded clinical data to filter the high dimensional data
tmpData = merge(allObservations$subjectInfo[c("subject.id","subject.inTrialId")],allObservations$observations,by="subject.id")
groups <- tmpData[c("subject.inTrialId","Subjects_Medical History_Tumor Characteristics_Tumor Subtypes")]
colnames(groups)<- c("subjects","tumorSubtypes")
groups <- groups[groups$"tumorSubtypes"!="None",]
View(groups)
expression_data <- expression_data[row.names(expression_data)%in% groups$subjects,]
View(expression_data)

## The following genes are a selection from all the genes used in the paper. Per subgroup a few genes have been used.
ERBB2 <- c("TLK2","MED24","MED1","ERBB2","GRB7") 
LumB <- c("ATP5G1","ERI3","CSDA","GGH","LAPTM4B","PRDX4","CCNE1")
LumA <- c("NAT1","SLC39A6","FOXA3","XBP1","GATA3","ESR1","PTP4A2","RERG","SCUBE2")
basal <- c("CXCL1","CDH3","ANXA8L1","KRT5","TRIM29","KRT17","CX3CL1","FZD8","CHI3L2","B3GNT5")
normal_breast <- c("PIK3R1","AKR1C1","ACSL1")
# gene list
genes <- c(ERBB2,LumA,LumB,basal,normal_breast)

# Transpose the expression_data to make it ready for generating a heatmap
mRNAdata <- t(expression_data[genes])

## expression data group labeling.
for (i in unique(groups$tumorSubtypes)){
  print(i)
  colnames(mRNAdata)[colnames(mRNAdata) %in% groups[groups$tumorSubtypes==i,"subjects"]] <- paste(
    i,colnames(mRNAdata)[colnames(mRNAdata) %in% groups[groups$tumorSubtypes==i,"subjects"]],sep="_")
}

## Generate a vector for the column corresponding to the groups.
colcolor<- colnames(mRNAdata)
colcolor[grep("Basal",colnames(mRNAdata))]<-"firebrick4"
colcolor[grep("Luminal A",colnames(mRNAdata))]<-"dodgerblue4"
colcolor[grep("Luminal B",colnames(mRNAdata))]<-"dodgerblue1"
colcolor[grep("ERBB2",colnames(mRNAdata))]<-"deeppink3"
colcolor[grep("Normal",colnames(mRNAdata))]<-"chartreuse4"

## generate heatmap ##
library(gplots)
heatmap.2(mRNAdata, Rowv = TRUE, Colv = TRUE, trace="none",
          col = greenred(800), ColSideColors = colcolor)

###############################################################################
# If you want to play more with the tranSMART Rclient possible options are:
# - Generate a larger heatmap
# - Try to remake the Kaplan Meier plot from the article
# - Play around with one of the many other interesting studies in the tranSMART instance
# For any questions please contact us at office@thehyve.nl or support@thehyve.nl
###############################################################################