library(dplyr)
library(tidyr)
2*pi
a <- 2*pi
a
b <- a <- 6
b
a
w <- c(1,2,9,8,7,5,pi)
w
v <- 1:10
v
colors <- c("red", "green", "blue", "white")
colors
v[4]
colors[4]
v > 6
letters
a <- letters[1:10]
a[v>6]
names(v) <- letters[1:10]
v
v["g"]
names(v) <- NULL
v[c(T,F)] #Addresses a question about subsetting with a conditional vector of a different length than v.
v[(c(T,F,F))] #When the length of v is not divisible by the length of the conditional vector used for subsetting.  
a[-c(2,4,6,8,10)]
m <- matrix(1:10, nrow=2, byrow=TRUE)
m
ma <- matrix(a, ncol=2)
ma
m[2,3]
ma[2,2]
m > 3
m[m > 3]
m[2,2:4]
m[2,]
m[,5]
m[,2:4]
m[1,]>3
m[,m[1,]>3]
Empl <- list(employee=c("Fred", "Barney"), spouse=c("Wilma", "Betty"), age=c(35,32), wage=c(10.5,8.7))
Empl
Empl$spouse
Empl$spouse[1]
Empl[[2]]
Empl[[2]][2]
Empl[3:4]
Empl[[3:4]]
Empl$employee <- c(Empl$employee, "Jetson")
Empl$employee
Empl <- c(Empl, new="Jetson")
unlist(Empl)
unlist(Empl, use.names=F)
fourcomp <- c(list(x=1:3, a=3:5), list(y=7:8, b=c(3,39)))
fourcomp
onecomp <- c(list(x=1:3, a=3:5), list(y=7:8, b=c(3,39)), recursive=T)
onecomp
exp <- factor(c("wt", "wt", "wt", "mut", "mut", "mut", "mut-tnf", "mut-tnf", "mut-tnf"))
exp
GeneExp <- c(2.1, 5.5, 14.2)
RespTNF <- c(F, T, T)
AffyId <- c("affyid1", "affyid2", "affyid3")
GeneName <- c("p53", "NFkB", "cMyc")
ExpData <- data.frame(GeneExp, RespTNF, GeneName, row.names=AffyId, stringsAsFactors=F)
ExpData
ExpData[1,]
ExpData[,2]
ExpData[,2] || ExpData[,2]
ExpData$GeneExp
row.names(ExpData)
ExpData[row.names(ExpData)=="affyid2",]
ExpDataTbl <- tbl_df(ExpData)
ExpDataTbl
?tbl_df
ExpDataTbl2 <- as.tbl(ExpData)
ExpDataTbl2
ExpDataRowNames <- tibble::rownames_to_column(ExpData, var = "affyId")
ExpDataRowNames
ExpDataTbl <- tbl_df(ExpDataRowNames)
ExpDataTbl
ExpDataTbl %>% filter(affyId == "affyid2")
ExpDataTbl %>% select(RespTNF)
ExpDataTbl %>% mutate(log2Exp = log2(GeneExp))
ExpDataTbl %>% mutate(log2Exp = log2(GeneExp)) %>% arrange(desc(GeneExp))
ExpDataTbl %>% summarize(mean(GeneExp))
ExpDataTbl %>% mutate(log2Exp = log2(GeneExp)) %>% summarize(mean(log2Exp))
typeof(a)
typeof(GeneExp)
attributes(ExpData)
ExpData[ExpData$RespTNF,]
mode(v)
length(v)
mode(ma)
length(ma)
mode(Empl)
length(Empl)
mode(ExpData)
length(ExpData)
mode(RespTNF)
length(RespTNF)
names(Empl)
names(ExpData)
names(ma)
dim(ma)
dim(ExpData)
dim(v)
dim(Empl)
nrow(ma)
ncol(ma)
help(t.test)
x <- rnorm(10)
y <- rnorm(10)
txy <- t.test(x, y , "greater")
names(txy)
txy$p.value
mode(txy)
t.test(x=x, alternative="greater", y=y)
m
tm <- t(m)
tm
write(tm, file="matrix.txt", ncol=5, sep="\t")
m_readin <- m <- matrix(scan("matrix.txt"), byrow=T, ncol=5)
m_readin
write.table(ExpData, file="ExpData.txt", quote=F, sep="\t", row.names=T, col.names=T) 
df <- read.table("ExpData.txt", header=T, as.is=T, sep="\t")
df
x <- seq(0,1,by=0.01)
y <- x + rnorm(length(x), mean=0, sd=0.1)
png(filename="linear_scatter.png")
plot(x, y, xlab="x", ylab="y", main="Linear Scatter", xlim=c(-.2,1), ylim=c(-.2,1), pch=18, col="red")
lines(x,x,col="blue")
graphics.off()
r <- rnorm(1000, mean=2, sd=1)
bmp(filename="hist.bmp")
hist(r, main="Hist of Rand Var r ~ N(2,1)")
graphics.off()
jpeg(filename="density.jpeg")
plot(density(r), main="Density of Rand Var r ~ N(2,1)")
graphics.off()
png("boxplot.png")
boxplot(r, main="Boxplot of r ~ N(2,1)")
graphics.off()
png("qqplot.png")
qqnorm(r, main="qqnorm of r ~ N(2,1)")
qqline(r)
graphics.off()
length(v)
if (length(v)>10) {long <- TRUE; variance <- var(v)} else {long <- FALSE; variance <- NA}
long
variance
square_root <- numeric()
for (i in 1:length(v)) {
  square_root <- c(square_root, sqrt(v[i]))
}
square_root
x<-1; i<-1; sum_sqrt<-0; while (square_root[i] <= 2) {sum_sqrt <- sum_sqrt + square_root[i]; i<-i+1}
sum_sqrt
med_mean <- function(x) {if (length(x) > 10) { mean(x) } else { median(x) }} 
data1 <- c(rnorm(6), 10)
med_mean(data1)
data1 <- c(rnorm(1000), 10) 
med_mean(data1) 
gt1 <- function(x) (sum(x>1))
data1 <- data.frame()
for (i in 1:100) {data1 <- rbind(data1, rnorm(10))}
dim(data1)
head(data1)
data1[1,]
rows <- apply(data1, 1, gt1)
columns <- apply(data1, 2, gt1)
length(rows)
length(columns)
names(columns) <- NULL; columns
highlow <- function(x) {if (x>10) {print("high")} else {print("low")}}
z <- 5:15
highlow(z) #This was associated with a question in class.  Doesn't work!  Scalar evaluation.
ifelse(z>10, "high", "low") #This works!  Element wise evaluation. 
