############
#
# collate_compounds() - R script for cpp3d
#
# AMI (rerun with new method, using AMI phyloseq data object from 2024 bioenergetic mapping study)
# Craig Liddicoat - Flinders University
# Running on Pawsey Setonix
############

# Add a new path
.libPaths(c("/software/projects/pawsey1216/cliddicoat/setonix/2024.05/r/4.4.1",
            "/software/projects/pawsey1216/cliddicoat/setonix/2024.05/r/4.3", .libPaths()))

R.Version()

# load packages
#library(readxl); packageVersion("readxl")
library(parallel); packageVersion("parallel")
library(doParallel); packageVersion("doParallel")
library(dplyr); packageVersion("dplyr")
library(stringr); packageVersion("stringr")
library(phyloseq); packageVersion("phyloseq") # '1.44.0'

no_forks <- 4

workdir <- "/scratch/pawsey1216/cliddicoat/ami"
setwd(workdir)
temp_dir <- "/scratch/pawsey1216/cliddicoat/ami/working"

this_study <- "-AMI"
phy <- readRDS("phy-phyloseq-object-AMI.RDS") # phyloseq object for SUPER-FOCUS functions
dat <- readRDS("dat.cpd-long-all-samps-cpp3d-ami.rds")

str(dat)

sum(dat$cpd_rel_abun_norm)
sum(dat$cpd_rel_abun_norm)/nsamples(phy) # average functional relative abundance per sample

length(which(is.na(dat$cpd_rel_abun_norm))==TRUE)
length(which( dat$cpd_rel_abun_norm > 0) == TRUE)
length(which( dat$cpd_rel_abun_norm == 0) == TRUE)

# so this step does collect some zero relative abundances from the 'otu-table' 

names(dat)
length(unique(dat$cpd_id)) # 


## Collate compounds within each sample 

unique_cpd <- unique(dat$cpd_id)
samp_names <- sample_names(phy)

collate_compounds <- function(dat.cpd, unique_cpd, samp) {
  #i<-1
  #samp = samp_names[i]
  #dat.cpd = dat[which(dat$sample == samp_names[i]), ]
  
  this_samp <- samp
  
  cpd_data <- data.frame(cpd_id = unique_cpd, sample=this_samp, cpd_rel_abun=NA)
  
  for (c in 1:length(unique_cpd)) {
    #c<-1
    this_cpd <- unique_cpd[c]
    sel.cpd <- which(dat.cpd$cpd_id == this_cpd)
    
    if (length(sel.cpd) >=1) {
      cpd_data$cpd_rel_abun[c] <- sum(dat.cpd$cpd_rel_abun_norm[sel.cpd])
    }
    
  } # END all compounds
  
  saveRDS(object = cpd_data, file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") ) # use readRDS()
  
} # END


#time.start <- Sys.time()
#cl<-makeCluster( detectCores()-2 )
#cl<-makeCluster( no_cores )

# this makes clusters on Unix-like system (may need to use other alternative for Windows)
cl<-makeForkCluster(nnodes = no_forks)      # no of nodes will depend on your HPC facility
registerDoParallel(cl)

foreach(i=1:length(sample_names(phy)), .packages=c('phyloseq')) %dopar%
  collate_compounds(dat.cpd = dat[which(dat$sample == samp_names[i]), ], unique_cpd = unique_cpd, samp = samp_names[i])

stopCluster(cl)
#time.finish <- Sys.time()


# output 1
i<-1
this_samp <- sample_names(phy)[i]
dat.cpd.collate <- readRDS ( file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") ) # use readRDS()
head(dat.cpd.collate)

for ( i in 2:length(sample_names(phy)) ) {
  #i<-1
  this_samp <- sample_names(phy)[i]
  temp <- readRDS ( file = paste0(temp_dir,"/cpd_data.collate-",this_samp,".rds") ) # use readRDS()
  
  dat.cpd.collate <- rbind(dat.cpd.collate, temp)
  
  print(paste0("completed ",i))
}


str(dat.cpd.collate)

sum(dat.cpd.collate$cpd_rel_abun)

sum(dat.cpd.collate$cpd_rel_abun)/length(unique(dat.cpd.collate$sample))

saveRDS(object = dat.cpd.collate, file = "dat.cpd.collate-all-samps-cpp3d-ami.rds" )

# END
