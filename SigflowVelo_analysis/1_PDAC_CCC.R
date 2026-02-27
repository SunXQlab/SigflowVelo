# Load the relevant packages
library(Seurat)
library(dplyr)
library(CellChat)
library(SummarizedExperiment)
library(zellkonverter)

### Set the ligand-receptor database. Here we will use the "Secreted signalling" database for cell-cell communication (let's look at ECM-receptor in the future )
CellChatDB <- CellChatDB.human # The othe roption is CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # Other options include ECM-Receptor and Cell-Cell Contact

### Burkhardt et al. (2021)
scPDAC <- readRDS("data/scPDAC_anno.rds")

scPDAC.data.input <- GetAssayData(scPDAC, assay = "RNA", slot = "data")

scPDAC.meta.data <- scPDAC@meta.data

#### First by merged cell types ####
scPDAC.identity <- data.frame(group = scPDAC.meta.data$celltype, row.names = row.names(scPDAC.meta.data))

Sensitive.use <- rownames(scPDAC.meta.data)[scPDAC.meta.data$group == "Sensitive"] # extract the cell names from disease data
Resist.use <- rownames(scPDAC.meta.data)[scPDAC.meta.data$group == "Resist"] # extract the cell names from disease data

Sensitive.data.input <- scPDAC.data.input[, Sensitive.use]
Resist.data.input <- scPDAC.data.input[, Resist.use]

Sensitive.meta <- scPDAC.meta.data[scPDAC.meta.data$group == "Sensitive", ]
Resist.meta <- scPDAC.meta.data[scPDAC.meta.data$group == "Resist", ]

scPDAC.Sensitive.identity <- data.frame(group = Sensitive.meta$celltype, row.names = row.names(Sensitive.meta))
colnames(scPDAC.Sensitive.identity)[1] <- "celltype"
scPDAC.Resist.identity <- data.frame(group = Resist.meta$celltype, row.names = row.names(Resist.meta))
colnames(scPDAC.Resist.identity)[1] <- "celltype"

# Create the cellchat objects
scPDAC.Sensitive.cc <-createCellChat(object = Sensitive.data.input, do.sparse = T, meta = scPDAC.Sensitive.identity, group.by = "celltype")
levels(scPDAC.Sensitive.cc@idents) # show factor levels of the cell labels

scPDAC.Resist.cc <-createCellChat(object = Resist.data.input, do.sparse = T, meta = scPDAC.Resist.identity, group.by = "celltype")
levels(scPDAC.Resist.cc@idents) # show factor levels of the cell labels

SensitiveGroupSize <- as.numeric(table(scPDAC.Sensitive.cc@idents)) # Get the number of cells in each group
ResistGroupSize <- as.numeric(table(scPDAC.Resist.cc@idents)) # Get the number of cells in each group

# Set the databases
scPDAC.Sensitive.cc@DB <- CellChatDB.use 
scPDAC.Resist.cc@DB <- CellChatDB.use 

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
### Control
scPDAC.Sensitive.cc <- subsetData(scPDAC.Sensitive.cc) # We subset the expression data of signalling genes to save on computational cost
scPDAC.Sensitive.cc <- identifyOverExpressedGenes(scPDAC.Sensitive.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
scPDAC.Sensitive.cc <- identifyOverExpressedInteractions(scPDAC.Sensitive.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
scPDAC.Sensitive.cc <- projectData(scPDAC.Sensitive.cc, PPI.human) # Other option includes PPI.human We're told that we may have to comment these out.

scPDAC.Sensitive.cc <- computeCommunProb(scPDAC.Sensitive.cc, raw.use = FALSE, population.size = FALSE)
scPDAC.Sensitive.cc <- computeCommunProbPathway(scPDAC.Sensitive.cc) # Calculate the probabilities at the signalling level
scPDAC.Sensitive.cc <- aggregateNet(scPDAC.Sensitive.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

### Stimulated
scPDAC.Resist.cc <- subsetData(scPDAC.Resist.cc) # We subset the expression data of signalling genes to save on computational cost
scPDAC.Resist.cc <- identifyOverExpressedGenes(scPDAC.Resist.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
scPDAC.Resist.cc <- identifyOverExpressedInteractions(scPDAC.Resist.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
scPDAC.Resist.cc <- projectData(scPDAC.Resist.cc, PPI.human) # Other option includes PPI.human We're told that we may have to comment these out.

scPDAC.Resist.cc <- computeCommunProb(scPDAC.Resist.cc, raw.use = FALSE, population.size = FALSE)
scPDAC.Resist.cc <- computeCommunProbPathway(scPDAC.Resist.cc) # Calculate the probabilities at the signalling level
scPDAC.Resist.cc <- aggregateNet(scPDAC.Resist.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

saveRDS(scPDAC.Resist.cc, file = "./CellChat/Output/scPDAC_Resist_cellchat.rds")
saveRDS(scPDAC.Sensitive.cc, file = "./CellChat/Output/scPDAC_Sensitive_cellchat.rds")

scPDAC.Sensitive.net <- subsetCommunication(scPDAC.Sensitive.cc)
scPDAC.Resist.net <- subsetCommunication(scPDAC.Resist.cc)

write.csv(scPDAC.Sensitive.net, file = "./CellChat/Output/scPDAC_celltype_communications_Sensitive.csv")
write.csv(scPDAC.Resist.net, file = "./CellChat/Output/scPDAC_celltype_communications_Resist.csv")

#### Now by the original clusters ####
# scPDAC.identity <- data.frame(group = scPDAC.meta.data$seurat_clusters, row.names = row.names(scPDAC.meta.data))

# Sensitive.use <- rownames(burkhardt.meta.data)[burkhardt.meta.data$Condition == "Ctrl"] # extract the cell names from disease data
# stim.use <- rownames(burkhardt.meta.data)[burkhardt.meta.data$Condition == "IFNg"] # extract the cell names from disease data

# control.data.input <- burkhardt.data.input[, control.use]
# stim.data.input <- burkhardt.data.input[, stim.use]

# control.meta <- burkhardt.meta.data[burkhardt.meta.data$Condition == "Ctrl", ]
# stim.meta <- burkhardt.meta.data[burkhardt.meta.data$Condition == "IFNg", ]

# burkhardt.control.identity <- data.frame(group = control.meta$leiden, row.names = row.names(control.meta))
# burkhardt.stim.identity <- data.frame(group = stim.meta$leiden, row.names = row.names(stim.meta))

# # Create the cellchat objects
# burkhardt.control.cc <-createCellChat(object = control.data.input, do.sparse = T, meta = burkhardt.control.identity, group.by = "group")
# levels(burkhardt.control.cc@idents) # show factor levels of the cell labels

# burkhardt.stim.cc <-createCellChat(object = stim.data.input, do.sparse = T, meta = burkhardt.stim.identity, group.by = "group")
# levels(burkhardt.stim.cc@idents) # show factor levels of the cell labels

# controlGroupSize <- as.numeric(table(burkhardt.control.cc@idents)) # Get the number of cells in each group
# stimGroupSize <- as.numeric(table(burkhardt.stim.cc@idents)) # Get the number of cells in each group

# # Set the databases
# burkhardt.control.cc@DB <- CellChatDB.use 
# burkhardt.stim.cc@DB <- CellChatDB.use 

# # We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
# ### Control
# burkhardt.control.cc <- subsetData(burkhardt.control.cc) # We subset the expression data of signalling genes to save on computational cost
# burkhardt.control.cc <- identifyOverExpressedGenes(burkhardt.control.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
# burkhardt.control.cc <- identifyOverExpressedInteractions(burkhardt.control.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
# burkhardt.control.cc <- projectData(burkhardt.control.cc, PPI.human) # Other option includes PPI.human We're told that we may have to comment these out.

# burkhardt.control.cc <- computeCommunProb(burkhardt.control.cc, raw.use = FALSE, population.size = FALSE)
# burkhardt.control.cc <- computeCommunProbPathway(burkhardt.control.cc) # Calculate the probabilities at the signalling level
# burkhardt.control.cc <- aggregateNet(burkhardt.control.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

# ### Stimulated
# burkhardt.stim.cc <- subsetData(burkhardt.stim.cc) # We subset the expression data of signalling genes to save on computational cost
# burkhardt.stim.cc <- identifyOverExpressedGenes(burkhardt.stim.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
# burkhardt.stim.cc <- identifyOverExpressedInteractions(burkhardt.stim.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
# burkhardt.stim.cc <- projectData(burkhardt.stim.cc, PPI.human) # Other option includes PPI.human We're told that we may have to comment these out.

# burkhardt.stim.cc <- computeCommunProb(burkhardt.stim.cc, raw.use = FALSE, population.size = FALSE)
# burkhardt.stim.cc <- computeCommunProbPathway(burkhardt.stim.cc) # Calculate the probabilities at the signalling level
# burkhardt.stim.cc <- aggregateNet(burkhardt.stim.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

# burkhardt.control.net <- subsetCommunication(burkhardt.control.cc)
# burkhardt.stim.net <- subsetCommunication(burkhardt.stim.cc)

# write.csv(burkhardt.control.net, file = "../复现/data/burkhardt21_leiden_communications_Ctrl.csv")
# write.csv(burkhardt.stim.net, file = "../复现/data/burkhardt21_leiden_communications_IFNg.csv")

