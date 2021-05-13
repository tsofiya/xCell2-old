# WIP: Samples are assumed to be seq. For final version- add array option.

# librarys
library(celldex)

# Directories
working.dir= "D:/сосиш з/xcell2/xCell/Dev_scripts/"

source(paste0(working.dir,'xCell_train_functions.R'))


# -- read reference data sets -- # Needs to be changed to user databaseref
data_sets = c('blueprint')
data_types = c('seq')

families=read.table(file.path(working.dir,'cells_families.txt'),sep="\t",header=TRUE,row.names=1, as.is=TRUE)
dependencies = read.types.dependencies(paste0(working.dir,'types_dependecies.txt'))

#read DB: to change to user data
ref <- BlueprintEncodeData()

gsc = create.signatures(ref,dependencies)

score.ref.signatures(ref, working.dir)

