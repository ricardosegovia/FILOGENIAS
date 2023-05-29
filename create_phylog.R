##Tutorial para crear una filogenia de un listado de especies

#devtools::install_github("jinyizju/V.PhyloMaker")
library(V.PhyloMaker)

library(vegan)
library(ape)
library(phytools)
library(phangorn)
library(taxize)
library(rgbif)


# Abrimos un listado de especies para hacer un arbol filogenetico
andes<-read.csv("jonispps.csv")
head(andes)
dim(andes)

## V.phylomaker necesita un data frame con tres columnas llasdas: "species", "genus", "family", "species.relative", "genus.relative" 

#Necesitamos separar el genero y conocer la familia.
#En este caso usamos la funcion de correccion taxonomica de rgbif
andes_temp <- name_backbone_checklist(andes$unique.spps.) #taxa biodata_corrected
head(as.data.frame(andes_temp))
#Es posible que algunos de estos nombres no sean aceptados
table(andes_temp$status)
#Seleccionaremos solo aquellos que son nombres aceptados
andes_temp2 <- andes_temp[andes_temp$status == "ACCEPTED",]
table(andes_temp2$status)
head(as.data.frame(andes_temp2))

#seleccionamos solo aquellas columnas que necesitamos
species_list <- as.data.frame(andes_temp2[,c("species", "genus", "family")])
head(species_list)

#creamos las columnas de relatives que necesita V.Phylomaker
species_list[, "species.relative"] <- NA
species_list[, "genus.relative"] <- NA
dim(species_list)
write.csv(species_list, "species.csv")

## Creamos el arbol
tree<-phylo.maker(sp.list =species_list, tree = GBOTB.extended, nodes = nodes.info.1, scenarios="S3")
#comprbamos el tipo de archivo creado
str(tree)
class(tree)
str(tree$scenario.3)
class(tree$scenario.3)

#guardamos el arbol
write.tree(tree$scenario.3, "temp_tree.tre")

## Ahora editamos el arbol
tree<-read.tree("temp_tree.tre")
str(tree)
newnames<-tree$tip.label
newnames <- gsub('_', ' ', newnames)

tree$tip.label<-newnames

#el arbol tiene 343 especies
dim(species_list)
length(tree$tip.label)
# hubo 15 especies que no fueron ubicadas
cuales <- species_list[!species_list$species%in%tree$tip.label,]
cuales

write.tree(tree,"tree.tre")

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

#ahora vamos a escoger solo a las especies de Liliopsida
node<-fastMRCA(tree,"Alstroemeria andina",
    "Mastigostyla cyrtophylla")
node
class(node)
liliop<-extract.clade(tree,node)
plotTree(liliop,ftype="i", node.numbers=T)

#subseleccion por numero de nodo

pooideae<-extract.clade(liliop,as.integer(53))
plotTree(pooideae,ftype="i", node.numbers=T)

