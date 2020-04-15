#BY Cláudia Constantino

#Among the genes previously selected by regularization with Cox model at the first 4 days 
#after injury is being selected to the ones with higher variance.

library(ggplot2)
library(ggpubr)

#compare genes with other temporal microarrays
summary(genes.selected.t0 %in% genes.selected.t1)
summary(genes.selected.t0 %in% genes.selected.t4)
summary(genes.selected.t1 %in% genes.selected.t4)


genes.selected <- c(genes.selected.t0,genes.selected.t1, genes.selected.t4)
#delete repeated genes
genes.selected <- unique(genes.selected)
length(genes.selected)

genes.data.subset <- patient.data[,1]
subset <- patient.data[,genes.selected]
genes.data.subset <- cbind(genes.data.subset,subset)

#genes.data.subset has the genes select for the 797 patient observations 

summary(genes.data.subset)

genes.data.subset <- as.data.frame(genes.data.subset)


#write.table(genes.data.subset, "genessubset.txt")


#eliminate microarray without information - Microarray ID: 23586509 and PATIENT ID: 16832704 
#correspond to index 391 in genes.data.subset dataframe
genes.data.subset <- genes.data.subset[-c(391), ] 
microarray.time <- as.data.frame(patient.data$SAMPLE_STUDYSTART_DAYS)
microarray.time <- microarray.time[-c(391),]



library(data.table)
DT <- data.table(genes.data.subset)

# this will get you your variance for each column
DT[, sapply(.SD, function(x) list(var=var(x)))]

# adding a `by` argument will give you the groupings
DT[, sapply(.SD, function(x) list(var=var(x))), by=PATIENT_ID]

# If you would like to round the values: 
DT[, sapply(.SD, function(x) list(var=round(var(x), 3))), by=PATIENT_ID]

# If we want to add names to the columns 
wide <- setnames(DT[, sapply(.SD, function(x) list(var=round(var(x), 3))), by=PATIENT_ID], c("ID", sapply(names(DT)[-1], paste0, c(".VAR"))))

summary(wide)

wide <- as.matrix(wide)
#a vector with the mean of each column
meanVar <- colMeans(wide[,2:84])
max <- max(meanVar)

wide.subset <- as.data.frame(wide)
wide.subset[colMeans(wide.subset) <= 0.05*max] <- NULL
#wide.subset[colMeans(wide.subset) > mean] <- NULL

#to remove a part of the string (remove '.VAR' from gene name)
colwide <- as.character(colnames(wide.subset[2:13]))
colwide <- gsub('.VAR', '', colwide)


#select in genes.data.subset just the selected genes with higher variance
genes.data.subset2 <- genes.data.subset[,colwide]
genes.data.subset2 <- cbind(genes.data.subset2,microarray.time, genes.data.subset$PATIENT_ID)

colnames(genes.data.subset2)[which(names(genes.data.subset2) == "microarray.time")] <- "MICROARRAY_TIME"
colnames(genes.data.subset2)[which(names(genes.data.subset2) == "genes.data.subset$PATIENT_ID")] <- "PATIENT_ID"


#LINE PLOT for the 12 genes

plot1 <- ggplot(genes.data.subset2, aes(x=genes.data.subset2$MICROARRAY_TIME, 
                               y=genes.data.subset2$NM_013450, 
                               group=genes.data.subset2$PATIENT_ID)) +
  geom_line(color="grey45", size=0.5)+
  geom_point(color="darkslategray3", size=0.7) +
  xlab("Time [days]") + ylab("Gene 'NM_013450'")



#Log transformation
genes.data.log <- log1p(genes.data.subset2[1:12])
genes.data.log <- cbind(genes.data.log, genes.data.subset2[13:14])

#line plot for the first gene with log-transformed values 
plotlog1 <- ggplot(genes.data.log, aes(x=genes.data.log$MICROARRAY_TIME, 
                               y=genes.data.log$NM_013450, 
                               group=genes.data.log$PATIENT_ID)) +
  geom_line(color="grey45", size=0.5)+
  geom_point(color="darkslategray3", size=0.7) +
  xlab("Time [days]") + ylab("Gene 'NM_013450' [log]")


figure1 <- ggarrange(plot1, plotlog1,
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2)
figure1

#Create 12 genes line plots in a loop
colNames <- names(genes.data.subset2)[1:12]
plot_list = list()
for(i in colNames){
  n <- noquote(i)
  plt <- ggplot(genes.data.subset2, aes_string(x=genes.data.subset2$MICROARRAY_TIME, 
                               y=genes.data.subset2[,i], 
                               group=genes.data.subset2$PATIENT_ID)) +
    geom_line(color="grey45", size=0.5)+
    geom_point(color="darkslategray3", size=0.3) +
    xlab("Time [days]") + ylab(n)
  print(plt)
  plot_list[[i]] = plt
  Sys.sleep(2)
}

figure2 <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 4)
figure2
#########



#Represent 
DT2 <- data.table(genes.data.subset2[,1:13])

# this will get you your variance for each column
DT2[, sapply(.SD, function(x) list(mean=mean(x)))]
# adding a `by` argument will give you the groupings
DT2[, sapply(.SD, function(x) list(mean=mean(x))), by=MICROARRAY_TIME]
# If you would like to round the values: 
DT2[, sapply(.SD, function(x) list(mean=round(mean(x), 3))), by=MICROARRAY_TIME]
# If we want to add names to the columns 
wide2 <- setnames(DT2[, sapply(.SD, function(x) list(mean=round(mean(x), 3))), by=MICROARRAY_TIME], c("time", sapply(names(DT2)[-13], paste0)))
summary(wide2)

#plot just the mean trajectory
ggplot(wide2, aes(x=wide2$time, y=wide2$NM_013450, 
                                group=1)) +
  geom_line(color="red4", size=0.5)+
  geom_point(color="red4", size=0.7) +
  xlab("Time [days]") + ylab("Gene 'NM_013450'")


#join CENTROID to plot1
p <- ggplot() +
  # plot gene trajectories for 168 patients
  geom_line(data=genes.data.subset2, aes(genes.data.subset2$MICROARRAY_TIME, 
                                         y=genes.data.subset2$NM_013450, 
                                         group=genes.data.subset2$PATIENT_ID),
              colour="grey45", size=0.5) +
  geom_point(data=genes.data.subset2, aes(x=genes.data.subset2$MICROARRAY_TIME, 
                                          y=genes.data.subset2$NM_013450, 
                                          group=genes.data.subset2$PATIENT_ID), 
             colour="darkslategray3", size=0.7) + 
  #plot gene mean trajectory
  geom_line(data=wide2, aes(x=wide2$time, y=wide2$NM_013450, 
                            group=1),
              colour="red4", size=0.7) +
  xlab("Time [days]") + ylab("Gene 'NM_013450'")


#12 genes line plots with CENTROID 
wide2 <- as.data.frame(wide2)

#function to move first column to the end
movetolast <- function(data, move) {
  data[c(setdiff(names(data), move), move)]
}
wide2 <- movetolast(wide2, c("time"))

colNames <- names(genes.data.subset2)[1:12]
plot_list2 = list()
for(i in colNames){
  n <- noquote(i)
  plt2 <- ggplot() +
    # plot gene trajectories for 168 patients
    geom_line(genes.data.subset2, mapping = aes(x=genes.data.subset2$MICROARRAY_TIME, 
                                             y=genes.data.subset2[,i], 
                                             group=genes.data.subset2$PATIENT_ID),
              colour="grey45", size=0.5) +
    geom_point(genes.data.subset2, mapping = aes(x=genes.data.subset2$MICROARRAY_TIME, 
                                              y=genes.data.subset2[,i], 
                                              group=genes.data.subset2$PATIENT_ID), 
               colour="darkslategray3", size=0.7) + 
    #plot gene mean trajectory
    geom_line(data=wide2, mapping = aes(x=wide2$time, y=wide2[,i], 
                              group=1),
              colour="red4", size=0.7) +
    xlab("Time [days]") + ylab(n)
  print(plt2)
  plot_list2[[i]] = plt2
  Sys.sleep(2)
}


# Save plots to PDF. Create pdf where each page is a separate plot
pdf("trajectories_mean.pdf")
for (i in 1:12) {
  print(plot_list2[[i]])
}
dev.off()




