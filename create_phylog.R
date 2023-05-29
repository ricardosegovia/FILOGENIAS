##creating an exploratory phylogeny
#devtools::install_github("jinyizju/V.PhyloMaker")

library(V.PhyloMaker)
library(vegan)
library(ape)
library(phytools)
library(phangorn)
library(taxize)
library(rgbif)


andes<-read.csv("/home/ricardo/Documents/Colaboraciones/johnny_riffo/tesis/jonispps.csv")
colnames(andes)
head(andes)
dim(andes)

## creating the file requiered by V.phylomaker "SPECIES","GENUS","FAMILY"
andes_temp <- name_backbone_checklist(andes$unique.spps.) #taxa biodata_corrected
head(as.data.frame(andes_temp))
table(andes_temp$status)
andes_temp2 <- andes_temp[andes_temp$status == "ACCEPTED",]
table(andes_temp2$status)
head(as.data.frame(andes_temp2))

data <- as.data.frame(andes_temp2[c("species", "genus", "family")])

species_list <- data
head(species_list)
species_list[, "species.relative"] <- NA
species_list[, "genus.relative"] <- NA

write.csv(species_list, "species.csv")

## Creating tree

tree<-phylo.maker(sp.list =species_list, tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S3")
str(tree)
write.tree(tree$scenario.3, "temp_tree.tre")

## editing tree - matching names with the species file
tree<-read.tree("temp_tree.tre")
str(tree)
newnames<-tree$tip.label
newnames <- gsub('_', ' ', newnames)

tree$tip.label<-newnames

write.tree(tree,"tree.tre")


length(tree$tip.label)

## plotting tree
plot(tree,no.margin=TRUE,edge.width=2)
plot(unroot(tree),type="unrooted",no.margin=TRUE,lab4ut="axial",
    edge.width=2)

plotTree(tree,type="fan",fsize=0.7,lwd=1,
    ftype="i")

#reducir el arbol - Eliminar Asteraceae y Fabaceae
table(species_list$family)
head(species_list)
fam2remov <- c("Asteraceae", "Fabaceae")

species2remov <- species_list[species_list$family%in%fam2remov,]$species


tree_corto <-drop.tip(tree,species2remov)
plotTree(tree_corto,type="fan",fsize=0.7,lwd=1,
    ftype="i")  

plotTree(tree_corto,type="fan",fsize=0.7,lwd=1,
    ftype="i")      


node<-fastMRCA(tree,"Alstroemeria andina",
    "Mastigostyla cyrtophylla")
class(node)
liliop<-extract.clade(tree,node)
plotTree(liliop,ftype="i", node.numbers=T)

#subseleccion por numero de nodo

pooideae<-extract.clade(liliop,as.integer(53))
plotTree(pooideae,ftype="i", node.numbers=T)

