source("C:\\Users\\JAL\\Desktop\\PhyloMeth\\Correlation\\CorrelationFunctions.R")
phy <- get_study_tree("ot_485", "tree1")
trait1 <- sim.char (phy,0.3, 1, model="BM")
trait2 <- sim.char (phy, 0.2, 1, model="speciational")
q <- matrix(c(-0.3, 0.3, 0.5, -0.5), 2, 2, byrow=TRUE)
discretetrait1 <- sim.char (phy, q, 1, model="discrete")
discretetrait2 <- sim.char (phy, q, 1, model="discrete")

discrete.data <- cbind(discretetrait1, discretetrait2)
rownames(discrete.data) <- rownames(discretetrait1)
continuous.data <- cbind(trait1, trait2)
rownames(continuous.data) <- rownames(trait1)

cleaned.continuous <- CleanData(phy, continuous.data)
cleaned.discrete <- CleanData(phy, discrete.data)
tree <- cleaned.discrete$phy
moss <- cleaned.discrete$data
moss[moss==2]<-0
fern <- cleaned.continuous$data

VisualizeData(tree, fern)
VisualizeData(tree, moss)
contrasts.answer <- RunContrasts(tree, fern)
save(list=ls(), file="CorrelationsResults.rda")
pagel94.answer <- RunPagel94(tree, moss)
save(list=ls(), file="CorrelationsResults.rda")

linregrmodel <- Runphylolm(tree, fern)
summary(linregrmodel)


