#####################################################################
########## R script for analysis of DiaOutline output files##########
#####################################################################

# Install package
install.packages("Momocs")

# Load package
library(Momocs)


##### Outline analysis - Species data
lf <- list.files("E:/DiaOutlineDataset/", pattern = "\\.txt$",full.names=TRUE)  
lf1<-lf_structure(lf,  split = "E:/DiaOutlineDataset/", trim.extension = FALSE)
lf2<-data.frame(substr(lf1$V2, 1, 5))
names(lf2)[1] <- "Type"

coordinates <- import_txt(lf) 
allDiatomOutlines<-Out(coordinates, fac = lf2)

panel(allDiatomOutlines, fac="Type", names=TRUE)
# stack(allDiatomOutlines)
stack(coo_center(allDiatomOutlines))
# calibrate_harmonicpower(allDiatomOutlines) #Estimates the number of harmonics required for the four Fourier methods implemented in Momocs
allDiatomOutlines.f <- efourier(allDiatomOutlines, nb.h=32) #Num. of harmonics set to 32

### PCA: Principal Component Analysis
allDiatomOutlines.p <- PCA(allDiatomOutlines.f)

tiff("E:/Supplement_10AllPCA.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 500)
plot(allDiatomOutlines.p, 1, chull.filled=TRUE,stars=TRUE, title="All data PCA") #Supplement 10 in the paper.
dev.off()

plot(allDiatomOutlines.p,"Type", title="All data PCA")
plot(allDiatomOutlines.p, 1, chull=TRUE, pos.shp = "full_axes", abbreviate.labelsgroups = TRUE, points=FALSE, labelspoints = TRUE)

### Linear Discriminant Analysis
allDiatomOutlines.l <- LDA(allDiatomOutlines.p,1)

tiff("E:/Fig8AllLDA.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 500)
plot(allDiatomOutlines.l, chull.filled=TRUE, stars=TRUE, title="All data LDA") #Fig 8 in the paper.
dev.off()

plot(allDiatomOutlines.l, title="All data LDA") #Another plot version
plot(allDiatomOutlines.l, 1, chull=TRUE, pos.shp = "full_axes", abbreviate.labelsgroups = TRUE, points=FALSE, labelspoints = TRUE) #Another plot version

##### MANOVA
m<-MANOVA(allDiatomOutlines.p, 'Type')
mpw<-MANOVA_PW(allDiatomOutlines.p, "Type")

sink("E:/Manova.txt")
lapply(mpw, print)
sink() 



####### Outline analysis - Genera data
lf <- list.files("E:/DiaOutlineDataset/", pattern = "\\.txt$",full.names=TRUE)  
lf1<-lf_structure(lf,  split = "E:/DiaOutlineDataset/", trim.extension = FALSE)
lf2<-data.frame(substr(lf1$V2, 1, 2))
names(lf2)[1] <- "Type"

coordinates <- import_txt(lf) 
allDiatomOutlines<-Out(coordinates, fac = lf2)

allDiatomOutlines.f <- efourier(allDiatomOutlines, nb.h=32) #Num. of harmonics set to 32

### PCA: Principal Component Analysis
allDiatomOutlines.p <- PCA(allDiatomOutlines.f)

tiff("E:/Supplement_9GeneraPCA.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 500)
plot(allDiatomOutlines.p, 1, chull.filled=TRUE,stars=TRUE, title="All data Genera PCA") #Supplement 9 in the paper.
dev.off()

plot(allDiatomOutlines.p,"Type", title="All data Genera PCA")
plot(allDiatomOutlines.p, 1, chull=TRUE, pos.shp = "full_axes", abbreviate.labelsgroups = TRUE, points=FALSE, labelspoints = TRUE)

### Linear Discriminant Analysis
allDiatomOutlines.l <- LDA(allDiatomOutlines.p,1)

tiff("E:/Fig3Genera.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 500)
plot(allDiatomOutlines.l, chull.filled=TRUE, stars=TRUE, title="Genera LDA") #Fig 3 in the paper.
dev.off()

plot(allDiatomOutlines.l, title="All data Genera LDA") #Another plot version
plot(allDiatomOutlines.l, 1, chull=TRUE, pos.shp = "full_axes", abbreviate.labelsgroups = TRUE, points=FALSE, labelspoints = TRUE) #Another plot version

##### MANOVA
m<-MANOVA(allDiatomOutlines.p, 'Type')
mpw<-MANOVA_PW(allDiatomOutlines.p, "Type")

sink("E:/GeneraManova.txt")
lapply(mpw, print)
sink()  



####### Outline analysis - Cymbella ONLY 
lf <- list.files("E:/DiaOutlineDataset/", pattern = "Cy(.*)txt$",full.names=TRUE)  
lf1<-lf_structure(lf,  split = "E:/DiaOutlineDataset/", trim.extension = FALSE)
lf2<-data.frame(substr(lf1$V2, 1, 5))
names(lf2)[1] <- "Type"

coordinates <- import_txt(lf) 
cyDiatomOutlines<-Out(coordinates, fac = lf2)

panel(cyDiatomOutlines, fac="Type", names=TRUE)
# stack(cyDiatomOutlines)
stack(coo_center(cyDiatomOutlines))
# calibrate_harmonicpower(cyDiatomOutlines) #Estimates the number of harmonics required for the four Fourier methods implemented in Momocs
cyDiatomOutlines.f <- efourier(cyDiatomOutlines, nb.h=32) #Num. of harmonics set to 32

### PCA: Principal Component Analysis
cyDiatomOutlines.p <- PCA(cyDiatomOutlines.f)

tiff("E:/Supplement_6CymPCA.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 500)
plot(cyDiatomOutlines.p, 1, chull.filled=TRUE,stars=TRUE, title="Cymbella PCA") #Supplement 6 in the paper.
dev.off()

plot(cyDiatomOutlines.p,"Type", title="Cymbella PCA")
plot(cyDiatomOutlines.p, 1, chull=TRUE, pos.shp = "full_axes", abbreviate.labelsgroups = TRUE, points=FALSE, labelspoints = TRUE)

### Linear Discriminant Analysis
cyDiatomOutlines.l <- LDA(cyDiatomOutlines.p,1)

tiff("E:/Fig6CymLDA.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 500)
plot(cyDiatomOutlines.l, chull.filled=TRUE, stars=TRUE, title="Cymbella LDA") #Fig 6 in the paper.
dev.off()

plot(cyDiatomOutlines.l, title="Cymbella LDA") #Another plot version
plot(cyDiatomOutlines.l, 1, chull=TRUE, pos.shp = "full_axes", abbreviate.labelsgroups = TRUE, points=FALSE, labelspoints = TRUE) #Another plot version



####### Outline analysis - Gomphonema ONLY 
lf <- list.files("E:/DiaOutlineDataset/", pattern = "Go(.*)txt$",full.names=TRUE)  
lf1<-lf_structure(lf,  split = "E:/DiaOutlineDataset/", trim.extension = FALSE)
lf2<-data.frame(substr(lf1$V2, 1, 5))
names(lf2)[1] <- "Type"

coordinates <- import_txt(lf) 
goDiatomOutlines<-Out(coordinates, fac = lf2)

panel(goDiatomOutlines, fac="Type", names=TRUE)
# stack(goDiatomOutlines)
stack(coo_center(goDiatomOutlines))
# calibrate_harmonicpower(goDiatomOutlines) #Estimates the number of harmonics required for the four Fourier methods implemented in Momocs
goDiatomOutlines.f <- efourier(goDiatomOutlines, nb.h=32) #Num. of harmonics set to 32

### PCA: Principal Component Analysis
goDiatomOutlines.p <- PCA(goDiatomOutlines.f)

tiff("E:/Supplement_5GomPCA.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 500)
plot(goDiatomOutlines.p, 1, chull.filled=TRUE,stars=TRUE, title="Gomphonema PCA") #Supplement 5 in the paper.
dev.off()

plot(goDiatomOutlines.p,"Type", title="Gomphonema PCA")
plot(goDiatomOutlines.p, 1, chull=TRUE, pos.shp = "full_axes", abbreviate.labelsgroups = TRUE, points=FALSE, labelspoints = TRUE)

### Linear Discriminant Analysis
goDiatomOutlines.l <- LDA(goDiatomOutlines.p,1)

tiff("E:/Fig4GomLDA.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 500)
plot(goDiatomOutlines.l, chull.filled=TRUE, stars=TRUE, title="Gomphonema LDA") #Fig 4 in the paper.
dev.off()

plot(goDiatomOutlines.l, title="Gomphonema LDA") #Another plot version
plot(goDiatomOutlines.l, 1, chull=TRUE, pos.shp = "full_axes", abbreviate.labelsgroups = TRUE, points=FALSE, labelspoints = TRUE) #Another plot version



####### Outline analysis - Gyrosigma ONLY 
lf <- list.files("E:/DiaOutlineDataset/", pattern = "Gy(.*)txt$",full.names=TRUE)  
lf1<-lf_structure(lf,  split = "E:/DiaOutlineDataset/", trim.extension = FALSE)
lf2<-data.frame(substr(lf1$V2, 1, 5))
names(lf2)[1] <- "Type"

coordinates <- import_txt(lf) 
gyDiatomOutlines<-Out(coordinates, fac = lf2)

panel(gyDiatomOutlines, fac="Type", names=TRUE)
# stack(gyDiatomOutlines)
stack(coo_center(gyDiatomOutlines))
# calibrate_harmonicpower(gyDiatomOutlines) #Estimates the number of harmonics required for the four Fourier methods implemented in Momocs
gyDiatomOutlines.f <- efourier(gyDiatomOutlines, nb.h=32) #Num. of harmonics set to 32

### PCA: Principal Component Analysis
gyDiatomOutlines.p <- PCA(gyDiatomOutlines.f)

tiff("E:/Supplement_8GyrPCA.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 500)
plot(gyDiatomOutlines.p, 1, chull.filled=TRUE,stars=TRUE, title="Gyrosigma PCA") #Supplement 8 in the paper.
dev.off()

plot(gyDiatomOutlines.p,"Type", title="Gyrosigma PCA")
plot(gyDiatomOutlines.p, 1, chull=TRUE, pos.shp = "full_axes", abbreviate.labelsgroups = TRUE, points=FALSE, labelspoints = TRUE)

### Linear Discriminant Analysis
gyDiatomOutlines.l <- LDA(gyDiatomOutlines.p,1)

tiff("E:/Fig7GyrLDA.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 500)
plot(gyDiatomOutlines.l, chull.filled=TRUE, stars=TRUE, title="Gyrosigma LDA") #Fig 7 in the paper.
dev.off()

plot(gyDiatomOutlines.l, title="Gyrosigma LDA") #Another plot version
plot(gyDiatomOutlines.l, 1, chull=TRUE, pos.shp = "full_axes", abbreviate.labelsgroups = TRUE, points=FALSE, labelspoints = TRUE) #Another plot version



####### Outline analysis - Luticola ONLY 
lf <- list.files("E:/DiaOutlineDataset/", pattern = "Lu(.*)txt$",full.names=TRUE)  
lf1<-lf_structure(lf,  split = "E:/DiaOutlineDataset/", trim.extension = FALSE)
lf2<-data.frame(substr(lf1$V2, 1, 5))
names(lf2)[1] <- "Type"

coordinates <- import_txt(lf) 
luDiatomOutlines<-Out(coordinates, fac = lf2)

panel(luDiatomOutlines, fac="Type", names=TRUE)
# stack(luDiatomOutlines)
stack(coo_center(luDiatomOutlines))
# calibrate_harmonicpower(luDiatomOutlines) #Estimates the number of harmonics required for the four Fourier methods implemented in Momocs
luDiatomOutlines.f <- efourier(luDiatomOutlines, nb.h=32) #Num. of harmonics set to 32

### PCA: Principal Component Analysis
luDiatomOutlines.p <- PCA(luDiatomOutlines.f)

tiff("E:/Supplement_7LutPCA.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 500)
plot(luDiatomOutlines.p, 1, chull.filled=TRUE,stars=TRUE, title="Luticola PCA") #Supplement 7 in the paper.
dev.off()

plot(luDiatomOutlines.p,"Type", title="Luticola PCA")
plot(luDiatomOutlines.p, 1, chull=TRUE, pos.shp = "full_axes", abbreviate.labelsgroups = TRUE, points=FALSE, labelspoints = TRUE)

### Linear Discriminant Analysis
luDiatomOutlines.l <- LDA(luDiatomOutlines.p,1)

tiff("E:/Fig5LutLDA.tiff", height = 12, width = 17, units = 'cm', compression = "lzw", res = 500)
plot(luDiatomOutlines.l, chull.filled=TRUE, stars=TRUE, title="Luticola LDA") #Fig 5 in the paper.
dev.off()

plot(luDiatomOutlines.l, title="Luticola LDA") #Another plot version
plot(luDiatomOutlines.l, 1, chull=TRUE, pos.shp = "full_axes", abbreviate.labelsgroups = TRUE, points=FALSE, labelspoints = TRUE) #Another plot version
