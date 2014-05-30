source("R_functions.R") # functions
library(lattice)
install.packages("gridExtra") # to use function 'grid.arrange' for plotting
require(gridExtra) 

####### DATA cleaning
####### Run this once and for all
####### Once this is done, directly go to ANALYSIS
dat.chem = read.csv("WNV EDIII RABBIT chem analysis.csv", check.names = F, colClasses = "character")
dat.hem = read.csv("WNV EDIII RABBIT hem analysis.csv", check.names = F, colClasses = "character")
colnames(dat.chem)
colnames(dat.hem)
all.equal(colnames(dat.chem)[1:7], colnames(dat.hem)[1:7])
which(dat.chem[,1:7]!=dat.hem[,1:7], arr.ind = T)
rbind(dat.chem[1,1:7], dat.hem[1,1:7])
all.equal(dat.chem[-1,1:7], dat.hem[-1,1:7])
tempdat = cbind(dat.hem[,1:7], dat.chem[,8:ncol(dat.chem)], dat.hem[,8:ncol(dat.hem)])
tempdat[tempdat=="NO DATA" | tempdat=="(-)"]="NA"
tempdat[tempdat=='5.5.'] = 5.5
#Trt1, Trt2, Trt3 are unique for each subject, so should be all the same for day 0, 28, 42 within the same subject. 
tempdat[3*1:48-2, 5:7] = tempdat[3*1:48-1, 5:7] 
colnames(tempdat)=gsub("\\\xb5", "micro", colnames(tempdat))
colnames(tempdat)=gsub("\\*", "x", colnames(tempdat))
colnames(tempdat)=gsub(" \\%", "\\%", colnames(tempdat))
colnames(tempdat)=gsub("\\%", " \\%", colnames(tempdat))
colnames(tempdat)= gsub("\\#", "", colnames(tempdat))
# remove range of values from the column names, e.g., " (5.0-8.0)"
colnames(tempdat)= gsub("[[:space:]]\\([0-9]+(\\.[0-9]+)?\\-[0-9]+(\\.[0-9]+)?\\)", "", colnames(tempdat))
# put units in parentheses
colnames(tempdat)= gsub("\\(", "", colnames(tempdat))
colnames(tempdat)= gsub("\\)", "", colnames(tempdat))
colnames(tempdat)[-(1:3)] = as.vector(sapply(colnames(tempdat)[-(1:3)], 
            function(x) {y=strsplit(x, " ")[[1]]; print(y[length(y)]); gsub(y[length(y)], paste("(", y[length(y)],")", sep = ""), x, fixed = T)}))

# write.csv(tempdat, file = "WNV EDIII RABBIT analysis.csv", row.names = F) 
# For variable 'ALB-PS', a value of '5.5.' was changed to '5.5' (repeated decimal point) 
# For variable 'TCHO-PS', a value of '(-)' was changed to 'NO DATA' as other missing values. 

############################# 
## Analysis
## Run from here
############################# 
# Read data:
mydat = read.csv("WNV EDIII RABBIT analysis.csv", check.names = F) 
colnames(mydat)
# col 8-20 (13 variables) for chem analysis, col 21-40 (20 variables) for hem analysis: 
# [1] "Experiment "                "Group"                      "Animal ID"                 
# [4] "Immunization time (days)"   "WNV DIII Vaccine (microg)"  "Adjuvant 1 MP7-OH (microg)"
# [7] "Adjuvant 2 CpG (microg)"    "Total Protein-PS (g/dl)"    "TBIL-PS (mg/dl)"           
# [10] "ALB-PS (g/dl)"              "TCHO-PS (mg/dl)"            "IP-PS (mg/dl)"             
# [13] "Ca-PS (mg/dl)"              "CRE-PS (mg/dl)"             "GGT-PS (U/l)"              
# [16] "GLU-PS (mg/dl)"             "BUN-PS (mg/dl)"             "ALT-PS (U/l)"              
# [19] "ALP-PS (U/l)"               "BUN-CRE (Ratio)"            "WBC (x1000/microL)"        
# [22] "NUET CTs (x1000/microL)"    "NEU (%)"                    "LYM CTs (x1000/microL)"    
# [25] "LYM (%)"                    "MONO CTs (x1000/microL)"    "MONO (%)"                  
# [28] "EOS CTs (x1000/microL)"     "EOS (%)"                    "BASO CTs (x1000/microL)"   
# [31] "BASO (%)"                   "RBC (x10^6/microL)"         "HGB (g/dL)"                
# [34] "HCT (%)"                    "MCV (fL)"                   "MCH (pg)"                  
# [37] "MCHC (g/dL)"                "RDW (%)"                    "PLT (x1000/microL)"        
# [40] "MPV (fL)"        
#### NOTE:
# columns 8-20: serum chemistry
# columns 21-40: CBC (hematology) 
####

# Combination of doses for each Group  
mydat[mydat[,4]==0,c(2,5:7)][4*1:12,]
grpinfo = mydat[match(1:11, mydat[,2]), c(2,5:7)]
#     Group WNV DIII Vaccine (microg) Adjuvant 1 MP7-OH (microg) Adjuvant 2 CpG (microg)
# 1       1                        45                        0.0                       0 #A00
# 13      2                        45                      150.0                       0 #AD0
# 25      3                        45                       42.6                       0 #AB0
# 37      4                        45                        0.0                      60 #A0B
# 49      5                        45                       21.3                      30 #AAA
# 61      6                        45                       42.6                      60 #ABB
# 85      7                        90                        0.0                       0 #B00
# 97      8                        90                       42.6                      60 #BBB
# 109     9                       135                        0.0                       0 #C00
# 121    10                       135                       42.6                      60 #CBB
# 133    11                       135                       64.0                      90 #CCC
paste(grpinfo[,2], grpinfo[,3], grpinfo[,4], sep = "+")
############

############ Create new data sets: 
############ Changes from day 0 to 28 and from day 0 to 42
rawall = mydat
colnames(rawall) = c("Experiment", "Group", "Subject", "Time", "Trt1", "Trt2", "Trt3", paste("x", seq(ncol(rawall)-7), sep = "."))
newlevels = c(1,3,2,4:11)
newlevels = c(1,5,3,4,6,2,7:11)
rawall$Group = factor(rawall$Group, levels = newlevels, labels = paste(grpinfo[,2], grpinfo[,3], grpinfo[, 4], sep = "+")[newlevels])
### take differenceces between two time points
change.0to28 = data.frame(rawall[rawall$Time==0, c(1:7)], rawall[rawall$Time==28, -(1:7)]-rawall[rawall$Time==0, -(1:7)], check.names = F)
change.0to42 = data.frame(rawall[rawall$Time==0, c(1:7)], rawall[rawall$Time==42, -(1:7)]-rawall[rawall$Time==0, -(1:7)], check.names = F)
change.28to42 = data.frame(rawall[rawall$Time==0, c(1:7)], rawall[rawall$Time==42, -(1:7)]-rawall[rawall$Time==28, -(1:7)], check.names = F)

#### save p-values
summary.pv = data.frame(
# between day 0 to 28, overall treatment effect
pv28.AllEff.Fstat = Fstat1(change.0to28)["Fstat",],
pv28.AllEff.numdf = Fstat1(change.0to28)["numdf",],
pv28.AllEff.dendf = Fstat1(change.0to28)["dendf",],
pv28.AllEff.pRAW = Fstat1(change.0to28)["pv",],
pv28.AllEff.pBF = pmin(Fstat1(change.0to28)["pv",]*33, 1), 
# between day 0 to 28, adjuvant effect
pv28.AdjEff.Fstat = Fstat2(change.0to28)["Fstat",],
pv28.AdjEff.numdf = Fstat2(change.0to28)["numdf",],
pv28.AdjEff.dendf = Fstat2(change.0to28)["dendf",],
pv28.AdjEff.pRAW = Fstat2(change.0to28)["pv",],
# between day 0 to 42, overall treatment effect
pv42.AllEff.Fstat = Fstat1(change.0to42)["Fstat",],
pv42.AllEff.numdf = Fstat1(change.0to42)["numdf",],
pv42.AllEff.dendf = Fstat1(change.0to42)["dendf",],
pv42.AllEff.pRAW = Fstat1(change.0to42)["pv",],
pv42.AllEff.pBF = pmin(Fstat1(change.0to42)["pv",]*33, 1), 
# between day 0 to 42, adjuvant effect
pv42.AdjEff.Fstat = Fstat2(change.0to42)["Fstat",],
pv42.AdjEff.numdf = Fstat2(change.0to42)["numdf",],
pv42.AdjEff.dendf = Fstat2(change.0to42)["dendf",],
pv42.AdjEff.pRAW = Fstat2(change.0to42)["pv",]
)
rownames(summary.pv) = names(mydat)[-(1:7)]
summary.pv
write.csv(summary.pv, "table_pvalues_final.csv")

#### plot fitted values 
which(summary.pv[,5]<0.05) # 5  9 10 12 13
which(summary.pv[,14]<0.05) #2  5  7 10 12 14
# some parameters for plotting
theme.novpadding <-
   list(layout.heights =
        list(top.padding = -.5,
 	    main.key.padding = 0,
 	    key.axis.padding = 0,
 	    axis.xlab.padding = 0,
 	    xlab.key.padding = 0,
 	    key.sub.padding = 0,
 	    bottom.padding = -.5),
        layout.widths =
        list(left.padding = 0,
 	    key.ylab.padding = 0,
 	    ylab.axis.padding = 0,
 	    axis.key.padding = 0,
 	    right.padding = 0))
for (K in c(5,14)) {
impvar = which(summary.pv[,K]<0.05) # important variables (variables with significant effects)
pdf(paste("plot_fitted_", if (K==5) "chem" else "CBC", ".pdf", sep = ""), width = 6.5, height = 3*length(impvar))
par(mar = c(0,0,0,0))
myplot = vector("list", length(impvar))
for (j in seq_along(impvar)) {
  i = 7+impvar[j]
  src = rep(c("obs", "fitted"), c(nrow(rawall), 9))
  base = mean(rawall[rawall$Time==0, i], na.rm = T)
  coef28 = coef(lm(change.0to28[,i]~factor(change.0to28$Trt1)-1))
  coef42 = coef(lm(change.0to42[,i]~factor(change.0to42$Trt1)-1))
  out = c(rawall[,i], rep(base,3), base + coef28, base + coef42)
  newdata = cbind(rbind(rawall[,c("Time", "Trt1", "Subject")], 
                  data.frame(Time = rep(c(0,28,42),each = 3), Trt1 = rep(unique(rawall$"Trt1"),3), Subject = rep(c("1","2","3"),3))),
                out)
  print(newdata)
  myplot[[j]] = xyplot(out~Time|factor(Trt1), groups = interaction(newdata$Subject, src), newdata, 
                  lty = rep(c(1,3),c(51,51)),
                  ylab = colnames(mydat)[i], 
                  type = "o", 
                  col = rep(c("red",grey(.5)),c(51,51)), 
                  distribute.type = F,
                  as.table = TRUE, 
                  aspect = 1, 
                  layout = c(3,1),
                  scales=list(x=list(at=c(0,28,42))),
                  xlab = colnames(mydat)[4],
                  par.settings = theme.novpadding)
                  
}
do.call(grid.arrange, c(myplot, ncol = 1, clip = T))
dev.off()
}
  