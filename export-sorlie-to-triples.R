# Copyright 2016 The Hyve BV.
#
# This file is part of tranSMART R Client: R package allowing access to
# tranSMART's data via its RESTful API.
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
require("transmartRClient")
connectToTransmart("http://transmart-demo.thehyve.net/transmart")

study <- "GSE4382_ETRIKS" 

#clinical data
observations <- getObservations(study, as.data.frame = F)
extractObservationTriple <- function(observation) { list(observation$subject$inTrialId, observation$label, observation$value)}
triples.list <- lapply(observations, extractObservationTriple)
triples.df <- as.data.frame(t(sapply(triples.list, rbind)))

#HD data?
lung_log_intensities <- getHighdimData(study.name = study, concept.match = "Lung", projection = "log_intensity")