# References for adjusted p-values for multiple testing: 
# Ge, Dudoit, and Speed (2003) 'resampling-based multiple testing for microarray data analysis'
# Dudoit, Shaffer and Boldrick (2003) 'Multiple hypothesis Testing in Microarray Experiments'

# vector of F-test statistics (one for each response variable)
Fstat1 = function(data) {
	sapply(data[,-(1:7)], function(y) {
#		m1 = lm(y~0, data)
		m2 = lm(y~0+(factor(data$Trt1)+factor(data$Trt2)+factor(data$Trt3))^2, data)
#		A = anova(m1,m2)
		f = summary(m2)$fstatistic # value, numdf, dendf
    	c(
    	Fstat = unname(f[1]), # same as A$"F"[2],  
      	pv = unname(pf(f[1],f[2],f[3], lower.tail = F)), #same as A$"Pr(>F)"[2]
      	numdf = unname(f[2]),
      	dendf = unname(f[3])
      	)
    })	
} #test: Fstat1(change.0to28)
lm(y~factor(data$Trt1))
Fstat2 = function(data) {
	sapply(data[,-(1:7)], function(y) {
		m1 = lm(y~0+factor(data$Trt1), data)
		m2 = lm(y~0+(factor(data$Trt1)+factor(data$Trt2)+factor(data$Trt3))^2, data)
		A = anova(m1,m2)
    	c(
    	Fstat = A$"F"[2],  
      	pv = A$"Pr(>F)"[2],
      	numdf = A$Df[2],
      	dendf = A$Res.Df[2]
      	)
    })	
} #test: Fstat2(change.0to28)
Fstat3 = function(data) {
	sapply(data[,-(1:7)], function(y) {
		m1 = lm(y~0+factor(data$Trt1), data)
		m2 = lm(y~0+(factor(data$Trt1)+factor(data$Trt2)+factor(data$Trt3))^3, data)
		A = anova(m1,m2)
    	c(
    	Fstat = A$"F"[2],  
      	pv = A$"Pr(>F)"[2],
      	numdf = A$Df[2],
      	dendf = A$Res.Df[2]
      	)
    })	
} #test: Fstat3(change.0to28)


# resample row indexes (either permute or bootstrap)
resample = function(x, method = c("permutation", "bootstrap"), resample.within = F) {
	method = match.arg(method)
	if (resample.within) { # resampling within levels of x
		ind = split(seq_along(x), x)
		ind.resample = lapply(ind, sample, replace = method == "bootstrap")
		unlist(ind.resample)[order(unlist(ind))]
	} else { # unrestricted resampling
		sample(seq_along(x), replace = method == "bootstrap")	
	}
} # test: resample(x=change.0to28$Trt1, resample.within = T)

# resampling raw p-values
resample.pv = function(data, B = 10, resample.method = c("permutation", "bootstrap"), resample.within = F) {
	resample.method = match.arg(resample.method)
	if (resample.within) Fstat = Fstat2 else Fstat = Fstat1
	# F-test statistics for the original data
    F.orig = Fstat(data)["Fstat",]
	m = ncol(data)-7 # number of response variables (number of genes, hypotheses)
    F.new = matrix(0,m,B)
	for (b in seq_len(B)) {	
		repeat{
			newind = resample(data$Trt1, method = resample.method, resample.within = resample.within)
			newdata = data
			newdata[,-(1:7)] = newdata[newind,-(1:7)]
			# compute F-test statistics
			F.new[,b] = try(Fstat(newdata)["Fstat",])
			if (!inherits(F.new, "try-error")) break
			# if (all(!is.na(F.new[, b]))) break
		}
	}	
	# resampling raw p-values
	p = rep(0,m)
	for (i in 1:m) p[i] = mean(F.new[i,]>=F.orig[i])
	p
} #test: resample.pv(change.0to28, B=1000, resample.method = "bootstrap", resample.within = F)

# compute B resampling test statistics and use the quick sort algorithm for each row (response variable, hypothsis, gene) to get the B raw p-values
# will be used in minP.sd.new
resample.pv.byquicksort = function(data, B = 10, resample.method = c("permutation", "bootstrap"), resample.within = F) { 
	resample.method = match.arg(resample.method)
	if (resample.within) Fstat = Fstat2 else Fstat = Fstat1
	m = ncol(data)-7 # number of response variables (number of genes, hypotheses)
    F.new = matrix(0,m,B) # resampling F-statistics
	for (b in seq_len(B)) {	
		repeat{
			newind = resample(data$Trt1, method = resample.method, resample.within = resample.within)
			newdata = data
			newdata[,-(1:7)] = newdata[newind,-(1:7)]
			# compute F-test statistics
			F.new[, b] = try(Fstat(newdata)["Fstat",])
			if (!inherits(F.new, "try-error")) break
			# if (all(!is.na(F.new[, b]))) break
		}
	}	
	# resampling raw p-values using the quick sort altorithm 
	t(apply(F.new, 1, function(x) rank(-x, ties.method = "max")/length(x)))
} # test: resample.pv.quicksort(change.0to28)

#### resampling adjusted p-values
# 1. Westfall and Young (1993) step-down max T procedure for FWER control
# use this if test statistics are identically distributed: in our case, v13:v33 have many missing values with differnent degrees of freedom from those for v1:v12
maxT.sd	= function(data, B = 10, resample.method = c("permutation", "bootstrap"), resample.within = F) { 
	resample.method = match.arg(resample.method)
	if (resample.within) Fstat = Fstat2 else Fstat = Fstat1
	# F-test statistics for the original data
    F.orig = Fstat(data)["Fstat",]#; print(F.orig)
    r = order(-F.orig)#; print(r)
	m = ncol(data)-7 # number of response variables
    u = matrix(0,m,B)
	for (b in seq_len(B)) {	
		repeat{
			newind = resample(data$Trt1, method = resample.method, resample.within = resample.within)
			newdata = data
			newdata[,-(1:7)] = newdata[newind,-(1:7)]
			# compute F-test statistics
			F.new = try(Fstat(newdata)["Fstat",])#; print(F.new)
			if (!inherits(F.new, "try-error")) break
			#if (all(!is.na(F.new))) break
		}
		# compute successive maxima of the statistics
		u[m,b] = F.new[r[m]]
		for (j in (m-1):1) u[j,b] = max(u[j+1,b], F.new[r[j]]) # print(u[,b])
	}
	# center and scale within each b
	# F.orig = scale(F.orig) 
	# u = scale(u)# NaN if variance=0 (constant)
	# resampling adjusted p-values
	p.adj = rep(0,m)
	for (j in 1:m) {p.adj[r[j]] = mean(u[j,]>=F.orig[r[j]])
		 if (is.na(mean(u[j,]>=F.orig[r[j]]))) {print(u[j,]); print(F.orig[r[j]]); print(mean(u[j,]>=F.orig[r[j]]))}
	}
	# enforce the monotonicity constraints
	for (j in 2:m) p.adj[r[j]] = max(p.adj[r[j]], p.adj[r[j-1]])
	p.adj
} #test: maxT.sd(change.0to28, B=10, resample.method = "permutation", resample.within = F)

# 2a. The traditional double permutation altorithm for Westfall and Young (1993) step-down min P procedure for FWER control: very slow
minP.sd	= function(data, B = 10, resample.method = c("permutation", "bootstrap"), resample.within = F) { 
	resample.method = match.arg(resample.method)
	# p-values for the original data
    #pv.orig = Fstat(data)[2,] # non-resampling p-values
    pv.orig = resample.pv(data, B=B, resample.method = resample.method, resample.within = resample.within) # resampling p-values
    r = order(pv.orig)
	m = ncol(data)-7 # number of response variables
    q = matrix(0,m,B)
	for (b in seq_len(B)) {	
		repeat{
			newind = resample(data$Trt1, method = resample.method, resample.within = resample.within)
			newdata = data
			newdata[,-(1:7)] = newdata[newind,-(1:7)]
			# compute p-values
			#pv.new = Fstat(newdata)[2,] #NaN if constant (variance=0) # non-resampling p-values
			pv.new = resample.pv(data, B=B, resample.method = resample.method, resample.within = resample.within) # resampling p-values
			if (all(!is.na(pv.new))) break
		}
		# compute successive minima of the p-values
		q[m,b] = pv.new[r[m]]
		for (j in (m-1):1) q[j,b] = min(q[j+1,b], pv.new[r[j]]) 
	}
	# center and scale within each b
	 # pv.orig = scale(pv.orig) 
	 # q = scale(q)# NaN if variance=0 (constant)
	# resampling adjusted p-values
	p.adj = rep(0,m)
	for (j in 1:m) p.adj[r[j]] = mean(q[j,]<=pv.orig[r[j]])
	# enforce the monotonicity constraints
	for (j in 2:m) p.adj[r[j]] = max(p.adj[r[j]], p.adj[r[j-1]])
	p.adj
} #test: minP.sd(change.0to28, B=5, resample.method = "permutation", resample.within = F)

# 2b. A new permutation algorithim for Westfall and Young (1993) step-down min P procedure for FWER control: same as minP.sd.new but faster
minP.sd.new	= function(data, B = 10, resample.method = c("permutation", "bootstrap"), resample.within = F) { 
	resample.method = match.arg(resample.method)
	# p-values for the original data
    #pv.orig = Fstat(data)[2,] # non-resampling p-values
    pv.orig.unsorted = resample.pv(data, B=B, resample.method = resample.method, resample.within = resample.within) # resampling p-values
    r = rank(pv.orig.unsorted) # needed later to go back to the original order of variables
    ord = order(pv.orig.unsorted)
    # ordered p-values
    pv.orig = pv.orig.unsorted[ord]
    # order data by p-values
    data[, -(1:7)] = data[, -(1:7)][, ord]
	m = ncol(data)-7 # number of response variables
    # compute B resampling test statistics and use the quick sort algorithm for each row (response variable, hypothsis, gene) to get the B raw p-values
	p = resample.pv.byquicksort(data, B = B, resample.method = c("permutation", "bootstrap"), resample.within = F)
	# initialize a matrix of successive minima
	q = matrix(1,m+1,B) # q[m+1,b] should be 1 for b=1,...,B
	# initialize a vector of adjusted p-values
	p.adj = rep(0,m)
	for (i in m:1) {
		# update the successive minima
		for (b in 1:B) q[i,b] = min(q[i+1,b], p[i,b])
		# adjusted p-value
		p.adj[i] = mean(q[i,]<=pv.orig[i])
		# delete row i of p
		if (i>=2) p = p[-i,,drop = F]
		# delete row i+i of q		
		if (i>=2) q = q[-(i+1),,drop = F]
	}
	# enforce monotonicity of p.adj
	for (j in 2:m) p.adj[j] = max(p.adj[j], p.adj[j-1])
	# go back to the original order of variables (note: x[order(x)][rank(x)] or sort(x)[rank(x)] equals x)
	data.frame(p.unadj = pv.orig.unsorted, p.adj = p.adj[r]) 
} #test: minP.sd.new(change.0to28, B=5, resample.method = "permutation", resample.within = F)
# set.seed(111); p1.b=minP.sd.new(change.0to28, B=10^4, resample.method = "bootstrap", resample.within = F)
   # p.unadj  p.adj
# 1   0.9445 1.0000
# 2   0.5034 0.9991
# 3   0.8964 1.0000
# 4   0.4088 0.9991
# 5   0.9024 1.0000
# 6   0.4818 0.9991
# 7   0.3035 0.9974
# 8   0.4823 0.9991
# 9   0.1575 0.9708
# 10  0.5454 0.9991
# 11  0.9293 1.0000
# 12  0.4738 0.9991
# 13  0.0866 0.8944
# 14  0.2419 0.9946
# 15  0.1708 0.9747
# 16  0.2434 0.9946
# 17  0.3602 0.9991
# 18  0.2839 0.9974
# 19  0.7741 1.0000
# 20  0.6628 0.9991
# 21  0.1090 0.9366
# 22  0.1186 0.9390
# 23  0.9719 1.0000
# 24  0.8443 1.0000
# 25  0.2793 0.9974
# 26  0.2880 0.9974
# 27  0.2438 0.9946
# 28  0.3786 0.9991
# 29  0.8291 1.0000
# 30  0.3665 0.9991
# 31  0.4052 0.9991
# 32  0.1152 0.9372
# 33  0.3926 0.9991
# set.seed(112); p2.b=minP.sd.new(change.0to42, B=10^4, resample.method = "bootstrap", resample.within = F)
   # p.unadj  p.adj
# 1   0.9688 0.9999
# 2   0.0113 0.2766
# 3   0.9134 0.9999
# 4   0.3277 0.9958
# 5   0.2616 0.9934
# 6   0.2801 0.9945
# 7   0.0442 0.6697
# 8   0.1693 0.9657
# 9   0.4032 0.9958
# 10  0.3820 0.9958
# 11  0.8911 0.9999
# 12  0.6826 0.9958
# 13  0.5022 0.9958
# 14  0.1736 0.9657
# 15  0.0299 0.5483
# 16  0.1007 0.8923
# 17  0.4201 0.9958
# 18  0.0375 0.6156
# 19  0.0078 0.2092
# 20  0.0483 0.6925
# 21  0.3484 0.9958
# 22  0.0794 0.8418
# 23  0.2836 0.9945
# 24  0.1657 0.9657
# 25  0.3389 0.9958
# 26  0.5741 0.9958
# 27  0.4451 0.9958
# 28  0.3021 0.9945
# 29  0.4734 0.9958
# 30  0.8886 0.9999
# 31  0.5848 0.9958
# 32  0.1022 0.8923
# 33  0.9346 0.9999
# set.seed(113); p3.b=minP.sd.new(change.0to28, B=10^2, resample.method = "bootstrap", resample.within = T)
# set.seed(114); p3.p=minP.sd.new(change.0to42, B=10^2, resample.method = "permutation", resample.within = T)


# > p1
 # [1] 1.000 0.993 1.000 0.993 0.880 0.993 0.975 0.993 0.748 0.993 1.000 0.993 0.755 0.977 0.931
# [16] 0.977 0.993 0.986 0.999 0.998 0.864 0.954 1.000 1.000 0.986 0.984 0.978 0.993 1.000 0.990
# [31] 0.993 0.750 0.993
# set.seed(111); y=resample.pv(change.0to28, B=1000)
 # [1] 0.951 0.367 0.930 0.319 0.080 0.418 0.185 0.378 0.052 0.415 0.939 0.489 0.052 0.192 0.123
# [16] 0.170 0.322 0.212 0.705 0.564 0.080 0.143 0.958 0.796 0.228 0.221 0.193 0.326 0.769 0.290
# [31] 0.325 0.063 0.308
# set.seed(112); p2=minP.sd.new(change.0to28, B=1000, resample.method = "permutation", resample.within = T)
# dim(change.0to28)