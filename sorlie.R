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

# For transmartRClient installation see https://github.com/thehyve/RInterface

# Example steps to authenticate with, connect to, and retrieve data from tranSMART
# load package
require("transmartRClient")

# Connect to the Hyve's demo transmart instance with the connectToTransmart function.
#   You can log in as a guest user with username "user" and password "user". 
connectToTransmart("http://transmart-demo.thehyve.net/transmart")

# retrieve a list of the available studies in the database:
studies <- getStudies()
print(studies)

# to access the studies programmatically use for example: 
# study<-studies$id[1]
# for the examples below we will use GSE4382, a study that ???
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
print(allObservations$observations[1:12,1:7]) #hint: if you are using RStudio you can also use the function "View" to see the data in a more user-friendly table, e.g.: View(allObservations$observations)
print(allObservations$subjectInfo[1:12,])
print(allObservations$conceptInfo[1:10,])

# retrieve information about the concepts that are present for this study
concepts <- getConcepts(study)
print(concepts)  

# You can use the concept names or the concept links to only retrieve data for a subset of the variables:

# retrieve observations for the first concept which conceptname contains "Tumor Subtypes"
observations <- getObservations(study, concept.match = "Tumor Subtypes", as.data.frame = T)
print(observations$observations)

# retrieve observations belonging to the first two concepts by using the api.links contained in the getConcepts-result
observations <- getObservations(study, concept.links = concepts$api.link.self.href[c(1,2)])
observations$observations[1:10,] 


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
dataDownloaded$data[1:10,1:10]

#select gene expression data
expression_data<-dataDownloaded$data[,-c(1:6)]
dim(expression_data)
rownames(expression_data)<-dataDownloaded$data$patientId

#Make a heatmap
#If the dimensions of the expression_data table are large, you may want to create a subset of the data first. 

# select only the cases and controls (excluding the patients for which the lung disease is not specified). Note: in the observation table the database IDs 
# are used to identify the patients and not the patient IDs that are used in the gene expression dataset
# TODO Define case and control properly
cases <- allObservations$observations$subject.id[allObservations$observations$'Subjects_Medical History_Tumor Characteristics_Tumor Subtypes' %in% c('Basal', 'ERBB2', 'Luminal A', 'Luminal B')]
controls <- setdiff(allObservations$observations$subject.id, cases)

# now we have the database IDs for the patients, but we need to get the patient IDs. These can be retrieved from the subjectInfo table: 
subjectInfo <- allObservations$subjectInfo
patientIDsCase    <- subjectInfo$subject.inTrialId[ subjectInfo$subject.id %in% cases ] 
patientIDsControl <- subjectInfo$subject.inTrialId[ subjectInfo$subject.id %in% controls] 

patientSet <- c(patientIDsCase,patientIDsControl)

patientSet<-patientSet[which(patientSet %in% rownames(expression_data))]
#make a subset of the data based on the selected patientSet and the probelist, and transpose the table so that the rows now contain probe names
#??? TODO Select probse
probeNames <- colnames(expression_data)[1:100]
subset<-t(expression_data[patientSet,probeNames]) 

#for ease of recognition: append "Case" and "Control" to the patient names
colnames(subset)[colnames(subset)%in% patientIDsCase] <- paste(colnames(subset)[colnames(subset)%in% patientIDsCase],"Case", sep="_" )
colnames(subset)[colnames(subset)%in% patientIDsControl] <- paste( colnames(subset)[colnames(subset)%in% patientIDsControl] , "Control",sep= "_")

dim(subset)
heatmap(as.matrix(subset), scale = "row")
