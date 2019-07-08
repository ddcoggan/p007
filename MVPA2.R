# p007 exp2 MVPA analysis pipeline
# subject info --------------------------------------------------------------------
ages <- c(28,15,35,33,24,24,18,19,23,25,21,26,20,23,20,41,22,25,45,24,25,23,25,24,34) # 11 male
# packages ---------------------------------------------------------------
packages <- c('ggplot2', 'reshape2', 'plyr', 'multcomp', 'broom', 'lsr', 'ez', 'plotrix', 'Hmisc',
              'nlme', 'grid', 'ppcor', 'psych','ggdendro','ggrepel','R.matlab')
lapply(packages, require, character.only = TRUE)

# configuration -----------------------------------------------------------
projDir <- '/Volumes/ddc/research/projects/p007/exp2'
categories = c('bottle','chair','face','house','shoe')
comparisons = c('within','between')
all_conds = c('Ibottle','Ichair','Iface','Ihouse','Ishoe','Sbottle','Schair','Sface','Shouse','Sshoe')
imageTypes = c('identical','similar', 'cross')
columns = data.frame(identical = c(1:5,11:15,21:25,31:35,41:45), similar = c(56:60,66:70,76:80,86:90,96:100), cross = c(51:55,61:65,71:75,81:85,91:95))
within_cols5 <- c(1,7,13,19,25)
between_cols5 <- c(6,11,12,16,17,18,21,22,23,24)
between_cols10 <- c(11,21,22,31:33,41:44,51:55,61:66,71:77,81:88,91:99)


# functions ---------------------------------------------------------------

sig_code <- function(x){
  y <- c()
  y[x > .05] = ''
  y[x <= .05] = '*'
  y[x <= .01] = '**'
  y[x <= .001] = '***'
  return(y)
}

roundDF <- function(x){
  is.num <- sapply(x, is.numeric)
  x[is.num] <- lapply(x[is.num], round, 5)
  return(x)
}

# GIST similarity matrix ----------------------------------

GISTFile <- file.path(projDir, 'scripts/stimuli/fSI.mat')
dataGIST <- readMat(GISTFile)

# set up blank matrix
GISTmat <- matrix(data=NA,nrow=10,ncol=10)

# start by getting GISTs for condition A
for (cA in 1:10){
  cA_name = all_conds[cA] # get condition name
  if (startsWith(cA_name,'I')){ itA <- 'i' } else { itA <- 's' } # ascertain whether this is identical or similar condition
  
  # if identical, we can just use the index of the category to get the GIST vectors
  if (itA == 'i'){
    gists_cA <- dataGIST$fSI[[1]][cA,,]
    
    # if similar, we need to find the stimuli, get the filenames, and use these to extract the correct GIST vectors.
  } else {
    gists_cA <- matrix(data = NA, nrow = 36, ncol = 4096)
    images <- Sys.glob(file.path(projDir, 'stimuli/similar/2-final/mvpa', categories[cA-5], "*.png"))
    imNames <- as.list(dataGIST$fSI[[4]][c(seq(1,16566,6))])
    for (i in 1:length(images)){
      imName <- basename(images[i])
      imName <- substr(imName,5,nchar(imName))
      imNo <- match(imName, imNames)
      imGIST <- dataGIST$fSI[[3]][imNo,]
      gists_cA[i,] <- imGIST
    }
  }
  
  # now get GISTs for condition B
  for (cB in 1:10){
    cB_name = all_conds[cB] # get condition name
    if (startsWith(cB_name,'I')){ itB <- 'i' } else { itB <- 's' } # ascertain whether this is identical or similar condition
    
    # if identical, we can just use the index of the category to get the GIST vectors
    if (itB == 'i'){
      gists_cB <- dataGIST$fSI[[1]][cB,,]
      
      # if similar, we need to find the stimuli, get the filenames, and use these to extract the correct GIST vectors.
    } else {
      gists_cB <- matrix(data = NA, nrow = 36, ncol = 4096)
      images <- Sys.glob(file.path(projDir, 'stimuli/similar/2-final/mvpa', categories[cB-5], "*.png"))
      imNames <- as.list(dataGIST$fSI[[4]][c(seq(1,16566,6))])
      for (i in 1:length(images)){
        imName <- basename(images[i])
        imName <- substr(imName,5,nchar(imName))
        imNo <- match(imName, imNames)
        imGIST <- dataGIST$fSI[[3]][imNo,]
        gists_cB[i,] <- imGIST
      }
    }
    
    # run GIST correlations, leaving one image out
    # this is different if we are looking within or between conditions
    gistCors <- c()
    if (cA == cB){
      # within
      for (i_ind in 1:36){
        imIndices <- 1:36
        g_ind <- gists_cA[i_ind,] # gist of left out image
        g_grp <- colMeans(gists_cA[imIndices[-i_ind],]) # gist of remaining group
        gistCors[i_ind] <- cor(g_ind,g_grp)
      }
    } else {
      # between
      for (i_ind in 1:36){
        g_ind <- gists_cA[i_ind,] # gist of left out image
        g_grp <- colMeans(gists_cB[,]) # gist of remaining group (all images, in this case)
        gistCors[i_ind] <- cor(g_ind,g_grp)
      }
    }
    GISTmat[cA,cB] <- mean(gistCors)
  }
}

# average GISTmat across diagonal
GISTmat <- (GISTmat+t(GISTmat))/2
for (r in 1:10){ for (c in 1:10){ if (r < c){ GISTmat[r,c] <- NA }}}

# remove diagonal for plot
GISTmat_plot = GISTmat
for (r in 1:10){ for (c in 1:10){ if (r == c){ GISTmat_plot[r,c] <- NA }}}

# plot
resdir <- file.path(projDir, 'results', 'models')
df = data.frame(cond = all_conds, t(GISTmat_plot))
colnames(df)[2:11] = all_conds
plotdf = melt(df, id = "cond")
plotdf$cond = factor(plotdf$cond, levels = all_conds)
minmax = c(round(min(plotdf$value, na.rm=T),digits=2),round(max(plotdf$value, na.rm = T),digits=2))
pdf(file = sprintf('%s/matrix_GIST.pdf',resdir),
    bg="transparent",
    width=4, 
    height=3)
print(ggplot(plotdf, aes(x=cond,y=ordered(variable, levels = rev(all_conds)),fill=value))+
        geom_tile(na.rm=F)+
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(),
              axis.text.x  = element_text(colour = "black", size=12, angle = 90, vjust = .3, hjust = 1),
              axis.text.y  = element_text(colour = "black", size=12, margin = margin(0,0,5,0)),
              axis.ticks = element_blank(),
              axis.title.x = element_text(colour = "black", size=12, margin = margin(10,0,0,0)), 
              axis.title.y = element_text(colour = "black", size=12, margin = margin(0,10,0,0)),
              legend.title = element_text(colour = "black", size=12, angle = 270)) +
        labs(x = "category    |     image  ", y = "   image    |    category") +
        scale_fill_continuous(low="red",high="yellow", na.value = "white", limits = minmax, breaks=minmax, oob = scales::squish,
                              guide = guide_colorbar(nbin=100, ticks=F, barwidth = .5, barheight = 10.0, title.position="right", title = "Pearson's r", title.hjust=.5))+
        scale_x_discrete(labels = rep(categories, 2)) +
        scale_y_discrete(labels = rep(rev(categories), 2)) +
        coord_fixed())
dev.off()

# compare withins and betweens for exp1 all conditions and exp2 'image' condition.

# within and between GIST cors for each cluster from experiment 1 manually calculated from other mvpa script and hard-coded
withins_exp1 <- fisherz(c(0.85, 0.92, 0.89, 0.93, 0.96, 0.95, 0.90, 0.88, 0.91, 0.82))
betweens_exp1 <- fisherz(c(-0.13555556, -0.19777778,  0.06333333, -0.09555556, -0.09888889, -0.09888889, -0.09555556,  0.06333333, -0.19777778, -0.13555556))
WB_exp1 <- withins_exp1-betweens_exp1
# get same for 'image' conditions in exp2
withins_exp2 <- c()
for (c in 1:5){ withins_exp2[c] <- GISTmat[c+5,c+5]}
betweens_exp2 <- c()
betweens_exp2[1] <- mean(GISTmat[7:10,6])
betweens_exp2[2] <- mean(c(GISTmat[8:10,7],GISTmat[7,6]))
betweens_exp2[3] <- mean(c(GISTmat[9:10,8],GISTmat[8,6:7]))
betweens_exp2[4] <- mean(c(GISTmat[10,9],GISTmat[9,6:8]))
betweens_exp2[5] <- mean(GISTmat[10,6:9])
WB_exp2 <- c(fisherz(withins_exp2)-fisherz(betweens_exp2))

WB_test <- t.test(WB_exp1, WB_exp2, paired = F)
W_test <- t.test(withins_exp1, withins_exp2, paired = F)


# semantic similarity matrix -----------------------------------------------
semMat=t(data.matrix(read.csv(file.path(projDir, 'data/wordnet_pat.csv'), header=F)))
for(x in 1:10){for(y in 1:10){ if (x<y){semMat[x,y]=NA}}} # remove top right half of matrix

# plot
resdir <- file.path(projDir, 'results', 'models')
df = data.frame(cond = all_conds, t(semMat))
colnames(df)[2:11] = all_conds
plotdf = melt(df, id = "cond")
plotdf$cond = factor(plotdf$cond, levels = all_conds)
minmax = c(round(min(plotdf$value, na.rm=T),digits=2),round(max(plotdf$value, na.rm = T),digits=2))
pdf(file = sprintf('%s/matrix_sem.pdf',resdir),
    bg="transparent",
    width=4.4, 
    height=3.3)
print(ggplot(plotdf, aes(x=cond,y=ordered(variable, levels = rev(all_conds)),fill=value))+
        geom_tile(na.rm=F)+
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(),
              axis.text.x  = element_text(colour = "black", size=12, angle = 90, vjust = .3, hjust = 1),
              axis.text.y  = element_text(colour = "black", size=12, margin = margin(0,0,5,0)),
              axis.ticks = element_blank(),
              axis.title.x = element_text(colour = "black", size=12, margin = margin(10,0,0,0)), 
              axis.title.y = element_text(colour = "black", size=12, margin = margin(0,10,0,0)),
              legend.title = element_text(colour = "black", size=12, angle = 270)) +
        labs(x = "category    |     image  ", y = "   image    |    category") +
        scale_fill_continuous(low="red",high="yellow", na.value = "white", limits = c(0,1), breaks=c(0,1), oob = scales::squish,
                              guide = guide_colorbar(nbin=100, ticks=F, barwidth = .5, barheight = 11.5, title.position="right", title = "path similarity", title.hjust=.5))+
        scale_x_discrete(labels = rep(categories, 2)) +
        scale_y_discrete(labels = rep(rev(categories), 2)) +
        coord_fixed())
dev.off()

# within v between
semVect <- as.vector(t(semMat))
iw <- semVect[c(1,12,23,34,45)]
ib <- semVect[c(11,21,22,31,32,33,41,42,43,44)]
sw <- semVect[c(56,67,78,89,100)]
sb <- semVect[c(66,76,77,86,87,88,96,97,98,99)]
xw <- semVect[c(51,62,73,84,95)]
xb <- semVect[c(52,53,54,55,61,63,64,65,71,72,74,75,81,82,83,85,91,92,93,94)]

plot_df <- data.frame(imageType = c('category','category','image','image','cat v im','cat v im'), 
                      comparison = c('within','between'),
                      mean = c(mean(iw), mean(ib), mean(sw), mean(sb), mean(xw), mean(xb)),
                      se = c(std.error(iw), std.error(ib), std.error(sw), std.error(sb), std.error(xw), std.error(xb)))
plot_df$imageType <- factor(plot_df$imageType, levels = c('category','image','cat v im'))
plot_df$comparison <- factor(plot_df$comparison, levels = c('within','between'))

# make bar plot
pdf(file=sprintf('%s/WB_semantic.pdf', resdir), 
    bg="transparent",
    width=5, 
    height=3)
print(ggplot(plot_df, aes(x=imageType, y=mean)) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line.x = element_line(colour = "black", size = .5), 
              axis.line.y = element_line(colour = "black", size = .5), 
              axis.ticks = element_line(colour = "black", size = .5),
              axis.text.x  = element_text(colour = "black", size=10), 
              axis.text.y  = element_text(colour = "black", size=10), 
              axis.title.x = element_blank(), 
              axis.title.y = element_text(colour = 'black', size=13, margin = margin(0,10,0,0)),
              legend.title = element_blank(),
              strip.text = element_text(colour = 'black', size = 13)) +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se, fill = comparison), width=.3, position=position_dodge(.6), size = .5) +
        geom_bar(aes(fill = comparison), stat = "identity", position = "dodge", width = .6, colour = "black") +
        geom_hline(yintercept = 0) +
        labs(y = "path similarity") +
        scale_fill_manual(values = c('white','grey66')) +
        scale_y_continuous(expand = c(0, 0), limits = c(0,1.1), breaks=seq(0,1,.2)))
dev.off()

# run t-tests and print to file
statsFile <-  file.path(resdir,'WB_semantic.txt')
sink(statsFile)
cat('### t-tests for within v between by image type in semantic model ###\n')

cat('\n## identical ##\n')
cat(capture.output(t.test(iw,ib,paired = F)),sep='\n')
cat('\n# Cohens D ##\n')
cat(capture.output(cohensD(iw,ib)),sep='\n')

cat('\n## similar ##\n')
cat(capture.output(t.test(sw,sb,paired = F)),sep='\n')
cat('\n# Cohens D ##\n')
cat(capture.output(cohensD(sw,sb)),sep='\n')

cat('\n## cross ##\n')
cat(capture.output(t.test(xw,xb,paired = F)),sep='\n')
cat('\n# Cohens D ##\n')
cat(capture.output(cohensD(xw,xb)),sep='\n')
sink()

# analysis within each ROI ----------------------------------------------------------------

for (dataType in c('within','between')){
  if (dataType == 'within'){
    data <- data.frame(read.table(file.path(projDir, 'data', 'matrices_within.csv'), header = F))
    data <- aggregate(x = data, by = list(data$V1,data$V2,data$V3), FUN = mean)
    data <- data[,c(1:3,8:107)]
  }
  if (dataType == 'between'){
    data <- data.frame(read.table(file.path(projDir, 'data', 'matrices.csv'), header = F))
  }
  
  colnames(data)[1:3] <- c('subject', 'maskProject', 'region')
  regions = levels(data$region)
  maskProjects = levels(data$maskProject)
  gap <- 3 # number of variable in 'data' before matrices begin
  DF <- length(levels(data$subject))-1
  
  for (maskProj in maskProjects[1:2]){
    
    dataMaskProj <- droplevels(data[data$maskProject == maskProj,])
    regions = levels(dataMaskProj$region)
    
    for (region in regions){
      
      print(region)
      
      data_region=dataMaskProj[dataMaskProj$region==region,]
      resdir=file.path(projDir, 'results/MVPA',dataType, maskProj, region)
      dir.create(resdir, recursive = T, showWarnings = F)
      
      # PLOT GROUP MEAN MATRICES
      mean_matrix_long = colMeans(data_region[,4:103])
      mean_matrix = matrix(data = mean_matrix_long, ncol = 10, nrow = 10)
      mean_matrix = t(mean_matrix)
      
      # make large matrix figure containing all conditions
      df = data.frame(cond = all_conds, t(mean_matrix))
      colnames(df)[2:11] = all_conds
      plotdf = melt(df, id = "cond")
      plotdf$cond = factor(plotdf$cond, levels = all_conds)
      minmax = c(round(min(plotdf$value, na.rm=T),digits=2),round(max(plotdf$value, na.rm = T),digits=2))
      pdf(file = sprintf('%s/matrix_all.pdf',resdir),
          bg="transparent",
          width=4, 
          height=3)
      print(ggplot(plotdf, aes(x=cond,y=ordered(variable, levels = rev(all_conds)),fill=value))+
              geom_tile(na.rm=F)+
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    axis.text.x  = element_text(colour = "black", size=12, angle = 90, vjust = .3, hjust = 1),
                    axis.text.y  = element_text(colour = "black", size=12, margin = margin(0,0,5,0)),
                    axis.ticks = element_blank(),
                    axis.title.x = element_text(colour = "black", size=12, margin = margin(10,0,0,0)), 
                    axis.title.y = element_text(colour = "black", size=12, margin = margin(0,10,0,0)),
                    legend.title = element_text(colour = "black", size=12, angle = 270)) +
              labs(x = "identical         similar", y = "similar         identical") +
              scale_fill_continuous(low="red",high="yellow", na.value = "white", limits = minmax, breaks=minmax, oob = scales::squish,
                                    guide = guide_colorbar(nbin=100, ticks=F, barwidth = .5, barheight = 10.0, title.position="right", title = "Pearson's R", title.hjust=.5))+
              scale_x_discrete(labels = rep(categories, 2)) +
              scale_y_discrete(labels = rep(rev(categories), 2)) +
              coord_fixed())
      dev.off()
      
      # make matrix figures for identical, similar and cross
      minmax = c(-.35,.53) # hard-coded range for the matrices in ventral stream
      for (cond in imageTypes){
        cols = columns[,colnames(columns) == cond]
        mean_matrix = matrix(data = mean_matrix_long[cols], ncol = 5, nrow = 5)
        mean_matrix = t(mean_matrix)
        
        # average across diagonal only for 'cross'
        if (cond == 'cross'){
          temp_mat = matrix(data=NA, nrow=5, ncol=5)
          for (r in 1:5){
            for (c in 1:5){
              if (r==c){
                temp_mat[r,c]=mean_matrix[r,c]}
              if(r>c){
                temp_mat[r,c] = mean(c(mean_matrix[r,c], mean_matrix[r,c]))}
            }
          }
          mean_matrix = temp_mat
        }
        
        # convert to vector and keep for RSA (unused as too few data points in between-only matrix)
        #assign(paste0('vector_', cond), c(t(mean_matrix))[between_cols5])
        
        # plot
        df = data.frame(cond = categories, t(mean_matrix))
        colnames(df)[2:6] = categories
        plotdf = melt(df, id = "cond")
        plotdf$cond = factor(plotdf$cond, levels = categories)
        #minmax = c(round(min(plotdf$value, na.rm=T),digits=2),round(max(plotdf$value, na.rm = T),digits=2))
        pdf(file = sprintf('%s/matrix_%s.pdf',resdir, cond),
            bg="transparent",
            width=3, 
            height=2)
        print(ggplot(plotdf, aes(x=cond,y=ordered(variable, levels = rev(categories)),fill=value))+
                geom_tile(na.rm=F)+
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(), 
                      panel.background = element_blank(),
                      axis.text.x  = element_text(colour = "black", size=10, angle = 90, vjust = .3, hjust = 1),
                      axis.text.y  = element_text(colour = "black", size=10, margin = margin(0,0,5,0)),
                      axis.ticks = element_blank(),
                      axis.title.x = element_blank(), 
                      axis.title.y = element_blank(), 
                      legend.title = element_text(colour = "black", size=10, angle = 270)) +
                scale_fill_continuous(low="red",high="yellow", na.value = "white", limits = minmax, breaks=minmax, oob = scales::squish,
                                      guide = guide_colorbar(nbin=100, ticks=F, barwidth = .5, barheight = 6.8, title.position="right", title = "Pearson's R", title.hjust=.5))+
                coord_fixed())
        dev.off()
      }
      
      # before any statistical analysis, z transform r values
      
      data_regionZ <- data_region
      data_regionZ[,4:103] <- fisherz(data_regionZ[,4:103])
      
      # REPRESENTATIONAL SIMILARITY ANALYSIS
      
      # # correlation between image types (group mean matrices)
      # GimtypeBetween <- matrix(data = c(fisherz(vector_identical), fisherz(vector_similar), fisherz(vector_cross)), ncol = 3)
      # intermatrixCors <- rcorr(GimtypeBetween)
      ## N of matrix clearly too low for this to be a useful analysis, therefore use individual matrices.
      
      # correlation between image types (individual matrices)
      combos <- combn(imageTypes,2)
      indRSA <- data.frame(subject = rep(levels(data$subject), each = dim(combos)[2]), imageTypeA = combos[1,], imageTypeB = combos[2,], r = NA)
      crossBetas <- c()
      catBetas <- c()
      GISTctrlSem = c()
      SemCtrlGIST = c()
      GISTctrlSem2 = c()
      SemCtrlGIST2 = c()
      GISTctrlSem3 = c()
      SemCtrlGIST3 = c()
      GISTnorm = c()
      Semnorm = c()
      GISTnorm2 = c()
      SemNorm2 = c()
      GISTnorm3 = c()
      SemNorm3 = c()
      for (s in levels(data$subject)){
        identical <- data_regionZ[data_regionZ$subject == s, columns$identical[between_cols5]+gap]
        similar <- data_regionZ[data_regionZ$subject == s, columns$similar[between_cols5]+gap]
        cross1 <- matrix(data = as.numeric(data_regionZ[data_regionZ$subject == s, columns$cross+gap]), nrow = 5, ncol = 5, byrow = T) # get 'cross' as matrix
        cross2 <- (cross1+t(cross1))/2 # average x diagonal
        cross3 <- c(cross2)
        cross <- cross3[between_cols5]
        indMats <- data.frame(identical = as.numeric(identical), similar = as.numeric(similar), cross = cross)
        for(contrast in 1:dim(combos)[2]){
          indRSA$r[indRSA$subject == s & indRSA$imageTypeA == combos[1,contrast] & indRSA$imageTypeB == combos[2,contrast]] <- cor(indMats[,colnames(indMats)==combos[1,contrast]], indMats[,colnames(indMats)==combos[2,contrast]])
        }
        
        # reviewer-requested analysis: multiple regression of category matrix with category v image matrix and categorical model
        identicalAll <- as.numeric(data_regionZ[data_regionZ$subject == s, columns$identical[c(within_cols5, between_cols5)]+gap])
        crossAll <- cross3[c(within_cols5, between_cols5)]
        catModel <- c(rep(1,5),rep(0,10))
        df <- data.frame(identical = identicalAll, cross = crossAll, category = catModel)
        model <- lm(identical ~ category + cross, data = df)
        crossBeta <- coef(summary(model))["cross","Estimate"]
        catBeta <- coef(summary(model))["category","Estimate"]
        crossBetas <- append(crossBetas, crossBeta)
        catBetas <- append(catBetas, catBeta)
        
        # run partial correlations
        GISTmatIdentical <- c(t(GISTmat[1:5,1:5]))[c(within_cols5, between_cols5)]
        GISTctrlSemTest <- pcor.test(identicalAll, GISTmatIdentical, catModel)
        GISTctrlSem = c(GISTctrlSem, GISTctrlSemTest$estimate)
        SemCtrlGISTtest <-  pcor.test(identicalAll, catModel, GISTmatIdentical)                            
        SemCtrlGIST <- c(SemCtrlGIST, SemCtrlGISTtest$estimate)
        
        # run normal correlations
        GISTtest <- cor.test(identicalAll, GISTmatIdentical)
        GISTnorm = c(GISTnorm, GISTtest$estimate)
        Semtest <-  cor.test(identicalAll, catModel)                            
        Semnorm <- c(Semnorm, Semtest$estimate)
        
        # run step-wise regression
        
        # IDENTICAL #
        # semantic after removing GIST
        model = lm(identicalAll ~ GISTmatIdentical)
        GISTnormTest = summary(model)
        GISTnorm2 <- c(GISTnorm2,GISTnormTest$coefficients[2])
        residual = residuals(model)
        SemCtrlGISTtest = summary(lm(residual ~ catModel))
        SemCtrlGIST2 <- c(SemCtrlGIST2, SemCtrlGISTtest$coefficients[2])
        
        # GIST after removing semantic
        model = lm(identicalAll ~ catModel)
        SemNormTest = summary(model)
        SemNorm2 <- c(SemNorm2, SemNormTest$coefficients[2])
        residual = residuals(model)
        GISTctrlSemTest = summary(lm(residual ~ GISTmatIdentical))
        GISTctrlSem2 <- c(GISTctrlSem2, GISTctrlSemTest$coefficients[2])
        
        # SIMILAR #
        similarAll <- as.numeric(data_regionZ[data_regionZ$subject == s, columns$similar[c(within_cols5, between_cols5)]+gap])
        GISTmatSimilar <- c(t(GISTmat[6:10,6:10]))[c(within_cols5, between_cols5)]
        
        # semantic after removing GIST
        model = lm(identicalAll ~ GISTmatSimilar)
        GISTnormTest = summary(model)
        GISTnorm3 <- c(GISTnorm3,GISTnormTest$coefficients[2])
        residual = residuals(model)
        SemCtrlGISTtest = summary(lm(residual ~ catModel))
        SemCtrlGIST3 <- c(SemCtrlGIST3, SemCtrlGISTtest$coefficients[2])
        
        # GIST after removing semantic
        model = lm(identicalAll ~ catModel)
        SemNormTest = summary(model)
        SemNorm3 <- c(SemNorm3, SemNormTest$coefficients[2])
        residual = residuals(model)
        GISTctrlSemTest = summary(lm(residual ~ GISTmatSimilar))
        GISTctrlSem3 <- c(GISTctrlSem3, GISTctrlSemTest$coefficients[2])
      }
      
      indRSA$z <- fisherz(indRSA$r)
      indRSA$contrast <- factor(paste(indRSA$imageTypeA, indRSA$imageTypeB, sep = '\n'))
      crossTest <- t.test(crossBetas, mu = 0)
      catTest <- t.test(catBetas, mu = 0)
      modTest <- t.test(catBetas,crossBetas, paired = T)
      regDF <- data.frame(model = c('cross', 'category'), mean = c(mean(crossBetas),mean(catBetas)), se = c(std.error(crossBetas),std.error(catBetas)), sig = c(crossTest$p.value, catTest$p.value))
      
      # run t-test against zero
      pw_df <- data.frame(imageTypeA = combos[1,], imageTypeB = combos[2,], mean = NA, se = NA, t = NA, df = NA, CI_lo = NA, CI_hi = NA, p = NA, p.sig = NA, p.cor = NA, p.cor.sig = NA, CohensD = NA)
      for(contrast in 1:dim(combos)[2]){
        values <- indRSA$z[indRSA$imageTypeA == combos[1,contrast] & indRSA$imageTypeB == combos[2,contrast]]
        test <- t.test(values, mu = 0)
        pw_df[contrast,3:13] <- c(mean(values), std.error(values), test$statistic, test$parameter, test$conf.int[1], test$conf.int[2], test$p.value, NA, NA, NA, cohensD(values, mu = 0))
      }
      pw_df$p.sig <- sig_code(pw_df$p)
      pw_df$p.cor <- p.adjust(pw_df$p, method = 'BH')
      pw_df$p.cor.sig <- sig_code(pw_df$p.cor)
      pw_df <- roundDF(pw_df)
      
      # print results to file
      corsFile <-  file.path(resdir,'MVPA_cors.txt')
      sink(corsFile)
      cat('### correlation between image types for between-category elements ###\n')
      cat(capture.output(pw_df),sep='\n')
      cat('\n### reviewer-requested multiple regression ###\n')
      cat(capture.output(regDF),sep='\n')
      cat('\n## test between models ##')
      cat(capture.output(modTest),sep='\n')
      cat('\n### reviewer-requested partial correlation ###\n')
      cat('\n## GIST controlling for W v B ##\n')
      cat(capture.output(t.test(GISTctrlSem, mu = 0)), sep='\n')
      cat('\n## W v B controlling for GIST ##\n')
      cat(capture.output(t.test(SemCtrlGIST, mu = 0)), sep='\n')
      cat('\n### Tim-requested partial correlation ###\n')
      cat('\n## GIST  ##\n')
      cat(capture.output(t.test(GISTnorm, mu = 0)), sep='\n')
      cat('\n## W v B ##\n')
      cat(capture.output(t.test(Semnorm, mu = 0)), sep='\n')
      cat('\n### Tim-requested step-wise regression ###\n')
      cat('\n### IDENTICAL ###\n')
      cat('\n## GIST ##\n')
      cat(capture.output(t.test(GISTnorm2, mu = 0)), sep='\n')
      cat('\n## W v B ##\n')
      cat(capture.output(t.test(SemNorm2, mu = 0)), sep='\n')
      cat('\n## GIST controlling for W v B ##\n')
      cat(capture.output(t.test(GISTctrlSem2, mu = 0)), sep='\n')
      cat('\n## W v B controlling for GIST ##\n')
      cat(capture.output(t.test(SemCtrlGIST2, mu = 0)), sep='\n')
      cat('\n### SIMILAR ###\n')
      cat('\n## GIST ##\n')
      cat(capture.output(t.test(GISTnorm3, mu = 0)), sep='\n')
      cat('\n## W v B ##\n')
      cat(capture.output(t.test(SemNorm3, mu = 0)), sep='\n')
      cat('\n## GIST controlling for W v B ##\n')
      cat(capture.output(t.test(GISTctrlSem3, mu = 0)), sep='\n')
      cat('\n## W v B controlling for GIST ##\n')
      cat(capture.output(t.test(SemCtrlGIST3, mu = 0)), sep='\n')
      sink()
      
      # make bar plot showing values for each contrast
      plot_df <- aggregate(data = indRSA, z ~ contrast, FUN = mean)
      temp <- aggregate(data = indRSA, z ~ contrast, FUN = std.error)
      plot_df$se <- temp$z
      plot_df$contrast <- factor(plot_df$contrast, levels = c('identical\nsimilar', 'identical\ncross', 'similar\ncross')) # force order of contrasts
      plot_df$sig <- pw_df$p.cor.sig[c(order(levels(plot_df$contrast), levels(indRSA$contrast)))] # put significance codes in (row orders pw_df and plot_df differ, so use indRSA order to correct)
      levels(plot_df$contrast) <- c('identical\nsimilar', 'identical\nident v sim', 'similar\nident v sim') # change labels of contrasts
      colnames(plot_df)[2] <- 'mean'
      pdf(file=(file.path(resdir, sprintf("barplot_corsXimtype.pdf", cond))), 
          bg="transparent",
          width=2, 
          height=3)
      print(ggplot(plot_df, aes(x=contrast, y=mean)) + 
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(), 
                    axis.line.x = element_line(colour = "black", size = .5), 
                    axis.line.y = element_line(colour = "black", size = .5), 
                    axis.ticks = element_line(colour = "black", size = .5),
                    axis.text.x  = element_text(colour = "black", size=10, angle = 90), 
                    axis.text.y  = element_text(colour = "black", size=10), 
                    axis.title.x = element_blank(), 
                    axis.title.y = element_text(size=13, margin = margin(0,10,0,0))) +
              geom_errorbar(aes(ymin=mean-se, ymax=mean+se), stat = "identity", width=.3, position=position_dodge(.7), size = .5) +
              geom_text(aes(label = sig, y = mean+se+.03), colour = "black", size = 5) +
              geom_bar(stat = "identity", position = "dodge", width = .7, colour = "black", fill = 'grey90') +
              labs(y = "correlation (z)") +
              coord_cartesian(ylim=c(-.1,1.4)) +
              scale_y_continuous(breaks=seq(-.8,2,.4), expand=c(0,0)))
      dev.off()
      
      # RSA between 10x10 GIST and neural matrices
      GISTvect <- na.omit(as.vector(t(GISTmat_plot)))
      NeurVect <- as.vector(colMeans(data_regionZ[,4:103])[between_cols10])
      statsFile <- file.path(resdir, 'GIST_v_neural_stats.txt')
      sink(statsFile)
      cat('### RSA between neural and GIST ###\n\n')
      cat(capture.output(cor.test(GISTvect,NeurVect)), sep ='\n')
      sink()
      
      # scatterplot
      plot_df = data.frame(matrixA = NeurVect, matrixB = GISTvect)
      pdf(file=file.path(resdir, "scatter_GISTvNeural.pdf"), 
          bg="transparent",
          width=5, 
          height=5)
      print(ggplot(plot_df, aes(x=matrixA, y=matrixB)) + 
              geom_smooth(method = "lm", alpha = .15, size = 1, fill = 'dodgerblue3', col='dodgerblue3') + 
              geom_point(size = 3) +
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(), 
                    axis.line.x = element_line(colour = "black", size = 1), 
                    axis.line.y = element_line(colour = "black", size = 1), 
                    axis.ticks = element_line(colour = "black", size = 1), 
                    axis.text.x  = element_text(colour = "black", size=18), 
                    axis.text.y  = element_text(colour = "black", size=18), 
                    axis.title.x = element_text(size=20, margin = margin(10,0,0,0)), 
                    axis.title.y = element_text(size=20, margin = margin(0,10,0,0))) +
              labs(x = 'fMRI response (z)' , y = 'image similarity (z)'))
      dev.off()
      
      # Tim-requested multiple regression
      SemVect <- c(c(semVect)[between_cols10])
      gistBetas = c()
      semBetas = c()
      for (s in levels(data_regionZ$subject)){
        indNeurVect <- as.numeric(data_regionZ[data_regionZ$subject == s, (between_cols10+gap)])
        df <- data.frame(neural = indNeurVect, GIST = GISTvect, semantic = SemVect)
        model <- lm(neural ~ GIST + semantic, data = df)
        gistBeta <- coef(summary(model))["GIST","Estimate"]
        semBeta <- coef(summary(model))["semantic","Estimate"]
        gistBetas <- append(gistBetas, gistBeta)
        semBetas <- append(semBetas, semBeta)
      }
      sink(corsFile, append = T)
      cat('\n### Tim-requested mutliple regression supermatrix ###\n')
      cat('\n## GIST betas against zero ##\n')
      cat(capture.output(t.test(gistBetas, mu = 0)), sep='\n')
      cat('\n## semantic betas against zero ##\n')
      cat(capture.output(t.test(semBetas, mu = 0)), sep='\n')
      cat('\n## GIST betas against semantic betas ##\n')
      cat(capture.output(t.test(gistBetas, semBetas, paired = T)), sep='\n')
      sink()
      
      # PATTERN RELIABILITY ANALYSIS
      data_WB = data.frame(subject = rep(levels(data$subject), each=30), imageType = rep(imageTypes, each=10), category = rep(categories, each=2), comparison = c('within', 'between'), z = NA)
      
      # state which columns in data correspond to which conditions (first letter is image type, second is w/b)
      iw <- c(1,12,23,34,45)+gap # 1 per category
      ib <- matrix(data = c(11,21,31,41,11,22,32,42,21,22,33,43,31,32,33,44,41,42,43,44)+gap, nrow = 5, ncol = 4, byrow = T) # 4 per category
      sw <- c(56,67,78,89,100)+gap # 1 per category
      sb <- matrix(data = c(66,76,86,96,66,77,87,97,76,77,88,98,86,87,88,98,96,97,98,99)+gap, nrow = 5, ncol = 4, byrow = T) # 4 per category
      xw <- c(51,62,73,84,95)+gap # 1 per category
      xb <- matrix(data = c(61,71,81,91,52,53,54,55,52,72,82,92,61,63,64,65,53,63,83,93,71,72,74,75,54,64,74,94,81,82,83,85,55,65,75,85,91,92,93,94)+gap, nrow = 5, ncol = 8, byrow = T) # 8 per category
      
      # input values into data_WB
      for (subj in levels(data$subject)){
        for (c in 1:length(categories)){
          data_WB$z[data_WB$subject == subj & data_WB$category == categories[c] & data_WB$imageType == 'identical' & data_WB$comparison == 'within'] <- data_regionZ[data_regionZ$subject == subj, iw[c]]
          data_WB$z[data_WB$subject == subj & data_WB$category == categories[c] & data_WB$imageType == 'identical' & data_WB$comparison == 'between'] <- mean(as.numeric(data_regionZ[data_regionZ$subject == subj, ib[c,]]))
          data_WB$z[data_WB$subject == subj & data_WB$category == categories[c] & data_WB$imageType == 'similar' & data_WB$comparison == 'within'] <- data_regionZ[data_regionZ$subject == subj, sw[c]]
          data_WB$z[data_WB$subject == subj & data_WB$category == categories[c] & data_WB$imageType == 'similar' & data_WB$comparison == 'between'] <- mean(as.numeric(data_regionZ[data_regionZ$subject == subj, sb[c,]]))
          data_WB$z[data_WB$subject == subj & data_WB$category == categories[c] & data_WB$imageType == 'cross' & data_WB$comparison == 'within'] <- data_regionZ[data_regionZ$subject == subj, xw[c]]
          data_WB$z[data_WB$subject == subj & data_WB$category == categories[c] & data_WB$imageType == 'cross' & data_WB$comparison == 'between'] <- mean(as.numeric(data_regionZ[data_regionZ$subject == subj, xb[c,]]))
        }
      }
      
      # run ANOVA and pairwise comparisons
      model <- ezANOVA(data=data_WB, dv=z, wid=subject, within=.(imageType, category, comparison), detailed=TRUE)
      data_WB$condition <- paste(data_WB$imageType, data_WB$category, data_WB$comparison, sep = '_')
      pw = pairwise.t.test(data_WB$z, data_WB$condition, paired = T, p.adjust.method = 'BH')
      pw_df <- data.frame(condition = paste(rep(imageTypes, each = 5), rep(categories,3), sep = '_'), mean = NA, se = NA, t = NA, df = NA, CI_lo = NA, CI_hi = NA, p = NA, p.sig = NA, p.cor = NA, p.cor.sig = NA, CohensD = NA)
      for (c in categories){
        for (i in imageTypes){
          W <- data_WB$z[data_WB$category == c & data_WB$imageType == i & data_WB$comparison == 'within']
          B <- data_WB$z[data_WB$category == c & data_WB$imageType == i & data_WB$comparison == 'between']
          test = t.test(W,B, paired = TRUE)
          pw_df[pw_df$condition == paste(i,c,sep='_'),2:12] = c(mean(W-B), std.error(W-B), test$statistic, 24, test$conf.int[1], test$conf.int[2], pw$p.value[paste(i,c,'within', sep = '_'), paste(i,c,'between', sep = '_')], NA, NA, NA, cohensD(W,B, method = 'paired'))
        }
      }
      pw_df$p.sig <- sig_code(pw_df$p)
      pw_df$p.cor <- p.adjust(pw_df$p, method = 'BH')
      pw_df$p.cor.sig <- sig_code(pw_df$p.cor)
      pw_df <- roundDF(pw_df)
      
      # print to file
      anovaFile <-  file.path(resdir,'MVPA_anova.txt')
      sink(anovaFile)
      cat('### ANOVA ### \n')
      cat(capture.output(model), sep='\n')
      cat('\n## pairwise contrasts across comparison, within image type and category ##\n')
      cat(capture.output(pw_df), sep='\n')
      sink()
      
      # arrange data for bar plot
      data_WB$imageType <- factor(data_WB$imageType, levels = imageTypes)
      data_WB$comparison <- factor(data_WB$comparison, levels = comparisons)
      plot_df <- aggregate(data = data_WB, z ~ subject + imageType + category, FUN = diff) # get difference in W, B
      temp <- aggregate(data = plot_df, z ~ category + imageType, FUN = mean) # get means
      plot_df <- aggregate(data = plot_df, z ~ category + imageType, FUN = std.error) # get ses
      colnames(plot_df)[3] = 'se'
      plot_df$mean <- -(temp$z) # put means into same df as ses while changing from B-W to W-B
      plot_df$sig <- pw_df$p.cor.sig # add sig. labels
      levels(plot_df$imageType) <- c('category','image','category v image')
      
      # make bar plot
      pdf(file=sprintf('%s/WB_imtypeXcategory.pdf', resdir), 
          bg="transparent",
          width=7, 
          height=3)
      print(ggplot(plot_df, aes(x=category, y=mean)) + 
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(), 
                    axis.line.x = element_line(colour = "black", size = .5), 
                    axis.line.y = element_line(colour = "black", size = .5), 
                    axis.ticks = element_line(colour = "black", size = .5),
                    axis.text.x  = element_text(colour = "black", size=10), 
                    axis.text.y  = element_text(colour = "black", size=10), 
                    axis.title.x = element_blank(), 
                    axis.title.y = element_text(colour = 'black', size=13, margin = margin(0,10,0,0)),
                    strip.text = element_text(colour = 'black', size = 13)) +
              facet_grid(. ~ imageType) +
              geom_errorbar(aes(ymin=mean, ymax=mean+se), width=.3, position=position_dodge(.9), size = .5) +
              geom_bar(stat = "identity", position = "dodge", width = .7, fill = "gray90", colour = "black") +
              geom_text(aes(label = sig, y = mean+se+.02), colour = "black", size = 5) +
              geom_hline(yintercept = 0) +
              labs(y = "within minus between (z)") +
              scale_y_continuous(expand = c(0, 0), limits = c(0,1.6), breaks=seq(0,2,.2)))
      dev.off()
      
      # now do permutation testing version
      iwAD <- c(1,3,6,10,15)
      ibAD <- matrix(data = c(2,4,7,11,2,5,8,12,4,5,9,13,7,8,9,14,11,12,13,14), nrow = 5, ncol = 4, byrow = T)
      swAD <- c(21,28,36,45,55)
      sbAD <- matrix(data = c(27,34,52,51,27,35,43,52,34,35,44,53,42:44,54,51:54), nrow = 5, ncol = 4, byrow = T)
      xwAD <- c(16,23,31,40,50)
      xbAD <- matrix(data = c(17:20,22,29,37,46,22,24:26,17,30,38,47,29,30,32,33,18,24,39,48,37:39,41,19,25,32,49,46:49,20,26,33,41), nrow = 5, ncol = 8, byrow = T)
      
      permutations <- array(data = NA, dim = c(25,3,5,1000))
      originals <- array(data = NA,dim = c(25,3,5))
      for (s in 1:25){
        temp_mat <- fisherz(as.numeric(data_regionZ[s,4:103]))
        temp_mat <- temp_mat[!is.na(temp_mat)]
        for (c in 1:5){
          originals[s,1,c] <- temp_mat[iwAD[c]] - mean(temp_mat[ibAD[c,]])
          originals[s,2,c] <- temp_mat[swAD[c]] - mean(temp_mat[sbAD[c,]])
          originals[s,3,c] <- temp_mat[xwAD[c]] - mean(temp_mat[xbAD[c,]])
          for (p in 1:1000){
            shuf_mat <- sample(temp_mat, size = length(temp_mat), replace = F)
            permutations[s,1,c,p] <- shuf_mat[iwAD[c]] - mean(shuf_mat[ibAD[c,]])
            permutations[s,2,c,p] <- shuf_mat[swAD[c]] - mean(shuf_mat[sbAD[c,]])
            permutations[s,3,c,p] <- shuf_mat[xwAD[c]] - mean(shuf_mat[xbAD[c,]])
          }
        }
      }
      permutations.collapse <- apply(permutations, c(2:4), mean)
      originals.collapse <- apply(originals, c(2,3), mean)
      dataPerm <- data.frame(imageType = rep(imageTypes, each = 5), category = categories, p = NA, p.sig = " ", p.cor = NA, p.cor.sig = " ")
      for (i in 1:3){
        for (c in 1:5){
          dataPerm$p[(i-1)*5+c] <- (length(permutations.collapse[i,c,permutations.collapse[i,c,]>originals.collapse[i,c]])+1)/1000
        }
      }
      dataPerm$p.sig <- sig_code(dataPerm$p)
      dataPerm$p.cor <- p.adjust(dataPerm$p, method = 'fdr')
      dataPerm$p.cor.sig <- sig_code(dataPerm$p.cor)
      
      # add results to ANOVA text file
      sink(anovaFile, append = T)
      cat('\n# permutation test for each cluster #\n')
      cat(capture.output(dataPerm), sep = '\n')
      sink() # unsink text file
      
      # run regression between neural and GIST within category across image type
      # in this instance, GIST refers to the 5 distributions of GIST similarities for each category across identical/similar
      imageData <- read.csv(file.path(projDir, 'stimuli/mvpa_means.csv'))
      gistVals <- imageData$r
      
      betas <- c()
      for (subj in levels(data$subject)){
        subjVals <- as.numeric(data[data$subject == subj & data$region == region, columns$cross[seq(1,25,6)]+gap])
        Ldata <- data.frame(neural = subjVals, gist = gistVals)
        Lmodel <- lm(neural ~ gist, data = Ldata)
        betas <- append(betas, coef(summary(Lmodel))["gist","Estimate"])
      }
      
      plot_df <- data.frame(mean = NA, se = NA, t = NA, df = DF, CI_lo = NA, CI_hi = NA, p = NA, p.sig = NA, cohensD = NA)
      plot_df$mean <- mean(betas)
      plot_df$se <- std.error(betas)
      model_test <- t.test(betas,mu = 0)
      plot_df$t <- model_test$statistic
      plot_df$CI_lo <- model_test$conf.int[1]
      plot_df$CI_hi <- model_test$conf.int[2]
      plot_df$p <- model_test$p.value
      plot_df$cohensD <- cohensD(betas, mu = 0)
      plot_df$p.sig <- sig_code(plot_df$p)
      
      statsFile <- file.path(resdir, 'GIST_v_neural_stats.txt')
      sink(statsFile, append = T)
      cat('\n### regression between neural and GIST within category across image type ###\n\n')
      cat(capture.output(plot_df), sep ='\n')
      sink()
      
    }
  }
}


# analysis across wang regions ------------------------------------------
resdir <- file.path(projDir, 'results','MVPA','between','Wang_2015','each_region')
dir.create(resdir, recursive = T, showWarnings = F)
wang_regions <- c('V1d','V1v','V2d','V2v','V3d','V3v','hV4','VO1','VO2','PHC1','PHC2','V3a','V3b','LO1','LO2','hMT')
streams <- c('posterior','ventral','lateral')

# set empty dataframe to be populated with withins and betweens
plot_df <- data.frame(region = rep(wang_regions, each = 3), stream = rep(streams, c(21,12,15)), imageType = imageTypes, mean = NA, se = NA)
pw_df <- data.frame(region = rep(wang_regions, each = 3), stream = rep(streams, c(21,12,15)), imageType = imageTypes, mean = NA, se = NA, t = NA, df = NA, CI_lo = NA, CI_hi = NA, p = NA, p.sig = NA, p.cor = NA, p.cor.sig = NA, CohensD = NA)

for (region in wang_regions){
  
  data_region=data[data$region==region,]
  data_regionZ <- data_region
  data_regionZ[,4:103] <- fisherz(data_regionZ[,4:103])
  
  withinColsIT <- list(c(1,12,23,34,45)+gap, c(56,67,78,89,100)+gap,c(51,62,73,84,95)+gap)
  betweenColsIT <- list(c(11,21,22,31,32,33,41,42,43,44)+gap, c(66,76,77,86,87,88,96,97,98,99)+gap, c(52,53,54,55,61,63,64,65,71,72,74,75,81,82,83,85,91,92,93,94)+gap)

  for (it in 1:3){
    
    withins <- c()
    betweens <- c()
    
    for (subj in levels(data$subject)){
      
      # store means of withins and betweens
      within <- mean(as.numeric(data_regionZ[data_regionZ$subject == subj, withinColsIT[[it]]]))
      between <- mean(as.numeric(data_regionZ[data_regionZ$subject == subj, betweenColsIT[[it]]]))
      withins <- append(withins, within)
      betweens <- append(betweens, between)
      
    }
    
    test <- t.test(withins, betweens, paired = T)
    pw_df[pw_df$region == region & pw_df$imageType == imageTypes[it],4:14] <- c(mean(withins-betweens), std.error(withins-betweens), test$statistic, 24, test$conf.int[1], test$conf.int[2], test$p.value, NA, NA, NA, cohensD(withins, betweens, method = 'paired'))
    
    # get W > B and place in plotting dataframe
    plot_df$mean[plot_df$region == region & plot_df$imageType == imageTypes[it]] <- mean(withins-betweens)
    plot_df$se[plot_df$region == region & plot_df$imageType == imageTypes[it]] <- std.error(withins-betweens)
    
  }
  pw_df$p.sig <- sig_code(pw_df$p)
  pw_df$p.cor <- p.adjust(pw_df$p, method = 'BH')
  pw_df$p.cor.sig <- sig_code(pw_df$p.cor)
  pw_df <- roundDF(pw_df)

  # put results in text file
  sink(file.path(resdir, 'stats.txt'))
  cat('### W v B for each region ###\n')
  cat(capture.output(pw_df), sep = '\n')
  sink()
}
  
plot_df$region <- factor(plot_df$region, levels = wang_regions)
plot_df$stream <- factor(plot_df$stream, levels = streams)
plot_df$imageType <- factor(plot_df$imageType, levels = imageTypes)
plot_df$imageType <- revalue(plot_df$imageType, c("cross"="identical v similar"))
plot_df$label <- pw_df$p.cor.sig

pdf(file=file.path(resdir, 'barplot_WB1.pdf'),
    bg="transparent",
    width=8, 
    height=6)
print(ggplot(plot_df, aes(x=region, y=mean)) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line.x = element_line(colour = "black", size = .5), 
              axis.line.y = element_line(colour = "black", size = .5), 
              axis.ticks = element_line(colour = "black", size = .5),
              axis.text.x  = element_text(colour = "black", size=10), 
              axis.text.y  = element_text(colour = "black", size=10), 
              axis.title.x = element_blank(),
              axis.title.y = element_text(colour = 'black', size=13, margin = margin(0,10,0,0)),
              strip.text = element_text(colour = 'black', size = 13)) +
        geom_errorbar(aes(ymin=mean, ymax=mean+se), stat = "identity", width=.3, position=position_dodge(.7), size = .5) +
        geom_bar(stat = "identity", position = "dodge", width = .7, colour = "black", fill = "grey90") +
        facet_grid(imageType~stream, scales="free", space="free_x", shrink=T) +
        geom_text(aes(label = label, y = mean+se+.02), colour = "black", size = 5) +
        geom_hline(yintercept = 0) +
        labs(y = 'within minus between (z)') +
        scale_y_continuous(limits = c(0,.8), breaks=seq(0,.6,.3), expand=c(0,0)))
dev.off()


# alternative plot, with region as a fill

plot_df$label[plot_df$label %in% c('**','***')] = '*' # very wide plot so dont distinguish different significance thresholds

pdf(file=file.path(resdir, 'barplot_WB2.pdf'),
    bg="transparent",
    width=14, 
    height=4)
print(ggplot(plot_df, aes(x=region, y=mean, fill = imageType)) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line.x = element_line(colour = "black", size = .5), 
              axis.line.y = element_line(colour = "black", size = .5), 
              axis.ticks = element_line(colour = "black", size = .5),
              axis.text.x  = element_text(colour = "black", size=13), 
              axis.text.y  = element_text(colour = "black", size=13), 
              axis.title.x = element_blank(),
              axis.title.y = element_text(colour = 'black', size=15, margin = margin(0,10,0,0)),
              legend.text = element_text(colour = 'black', size = 13),
              legend.title = element_blank(),
              strip.text = element_text(colour = 'black', size = 13)) +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), stat = "identity", width=.3, position=position_dodge(.7), size = .5) +
        geom_bar(stat = "identity", position = "dodge", width = .7, colour = "black") +
        facet_grid(.~stream, scales="free", space="free_x", shrink=T) +
        geom_text(aes(label = label, x = region, y = mean+se+.02), position = position_dodge(width = .7), colour = "black", size = 5) +
        geom_hline(yintercept = 0) +
        scale_fill_manual(values = c('white','grey66','grey33')) +
        labs(y = 'within minus between (z)'))
#        scale_y_continuous(breaks=seq(ybreaks_min[model], ybreaks_max[model], ybreaks_step[model]), expand=c(0,0)))
dev.off()


# region matrix
Nregions = length(wang_regions)
all_mean_mats = matrix(data = NA, nrow = Nregions, ncol = 45) # set up blank matrix
rownames(all_mean_mats) = wang_regions
for (r in 1:Nregions){
  all_mean_mats[r,] <- fisherz(colMeans(data[data$region == wang_regions[r],between_cols10+3]))
}
inv_cors_mat <- data.frame(region = wang_regions, 1-(cor(t(all_mean_mats)))) # matrix of 1-R
inv_cors_mat[inv_cors_mat == 0] = NA
temp_plot_data <- melt(inv_cors_mat, id = "region")
temp_plot_data$region <- factor(temp_plot_data$region, levels = wang_regions) # change order of regions
temp_plot_data$variable <- factor(temp_plot_data$variable, levels = rev(wang_regions))
minmax=c(min(temp_plot_data$value, na.rm=T)-.01, max(temp_plot_data$value, na.rm=T)+.01)
rounded_minmax <- round(minmax, digits=2)
pdf(file = file.path(resdir, "region_matrix.pdf"),
    bg="transparent",
    width=7, 
    height=6)
print(ggplot(temp_plot_data, aes(x=region,y=variable,fill=value))+
        geom_tile(na.rm=F)+
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(),
              axis.text.x  = element_text(colour = "black", size=12, angle = 90, hjust = 1, vjust = .5),
              axis.text.y  = element_text(colour = "black", size=12),
              axis.title.x = element_blank(), 
              axis.title.y = element_blank(),
              legend.title = element_text(colour = "black", size=12, angle = 270),
              legend.text = element_text(colour = 'black', size = 12)) +
        scale_fill_continuous(low="yellow",high="red", na.value = "grey", limits = c(rounded_minmax), breaks=c(rounded_minmax),
                              guide = guide_colorbar(nbin=100, ticks=F, barwidth = 1, barheight = 26, title.position="right", title = "distance (1-r)", title.hjust=.5))+
        coord_fixed())
dev.off()

# hierarchical clustering
dist_mat <- dist(all_mean_mats, method = 'maximum') # matrix of euclidean distances
hc_data <- hclust(dist_mat, method = 'complete')
pdf(file = file.path(resdir, "hier_cluster.pdf"),
    bg="transparent",
    width=4, 
    height=7)
print(ggdendrogram(hc_data, rotate = T) +
        theme(axis.text.x = element_text(size = 12, colour = 'black', hjust = .5),
              axis.text.y = element_text(size = 12, colour = 'black'),
              axis.line.x = element_line(size = .5),
              axis.ticks.x = element_line(size = .5),
              axis.title.x = element_text(size = 15, colour = 'black')) +
        labs(y = 'maximum distance'))
dev.off()

# multi-dimensional scaling
fit = cmdscale(dist_mat, eig = T, k = 2)
x = fit$points[,1]
y = fit$points[,2]
mds_plot_data = data.frame(x = x, y = y)
pdf(file = file.path(resdir, "MDS.pdf"),
    bg="transparent",
    width=7, 
    height=7)
print(ggplot(mds_plot_data, aes(x=x,y=y))+
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", size = 2, fill = "transparent"),
              axis.text.x  = element_blank(),
              axis.text.y  = element_blank(),
              axis.ticks.x  = element_blank(),
              axis.ticks.y  = element_blank(),
              axis.title.x = element_blank(), 
              axis.title.y = element_blank(),
              legend.position="none") +
        geom_point(na.rm=F)+
        geom_text_repel(aes(label = rownames(mds_plot_data), size = 5),box.padding = unit(0.25, "lines")) +
        coord_fixed())
dev.off()

# region matrices, hierarchical clustering and MDS separately for each image type
for (it in c('identical', 'similar')){
  cols2use <- columns[between_cols5,colnames(columns) == it]
  all_mean_mats_it = matrix(data = NA, nrow = Nregions, ncol = 10) # set up blank matrix
  rownames(all_mean_mats_it) = wang_regions
  for (r in 1:Nregions){
    all_mean_mats_it[r,] <- fisherz(colMeans(data[data$region == wang_regions[r],cols2use+3]))
  }
  inv_cors_mat <- data.frame(region = wang_regions, 1-(cor(t(all_mean_mats_it)))) # matrix of 1-R
  inv_cors_mat[inv_cors_mat == 0] = NA
  temp_plot_data <- melt(inv_cors_mat, id = "region")
  temp_plot_data$region <- factor(temp_plot_data$region, levels = wang_regions) # change order of regions
  temp_plot_data$variable <- factor(temp_plot_data$variable, levels = rev(wang_regions))
  minmax=c(min(temp_plot_data$value, na.rm=T)-.01, max(temp_plot_data$value, na.rm=T)+.01)
  rounded_minmax <- round(minmax, digits=2)
  pdf(file = file.path(resdir, sprintf("region_matrix_%s.pdf", it)),
      bg="transparent",
      width=7, 
      height=6)
  print(ggplot(temp_plot_data, aes(x=region,y=variable,fill=value))+
          geom_tile(na.rm=F)+
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.background = element_blank(),
                axis.text.x  = element_text(colour = "black", size=12, angle = 90, hjust = 1, vjust = .5),
                axis.text.y  = element_text(colour = "black", size=12),
                axis.title.x = element_blank(), 
                axis.title.y = element_blank(),
                legend.title = element_text(colour = "black", size=12, angle = 270),
                legend.text = element_text(colour = 'black', size = 12)) +
          scale_fill_continuous(low="yellow",high="red", na.value = "grey", limits = c(rounded_minmax), breaks=c(rounded_minmax),
                                guide = guide_colorbar(nbin=100, ticks=F, barwidth = 1, barheight = 26, title.position="right", title = "distance (1-r)", title.hjust=.5))+
          coord_fixed())
  dev.off()
  
  # hierarchical clustering
  dist_mat <- dist(all_mean_mats_it, method = 'euclidean') # matrix of euclidean distances
  hc_data <- hclust(dist_mat, method = 'complete')
  pdf(file = file.path(resdir, sprintf("hier_cluster_%s.pdf", it)),
      bg="transparent",
      width=4, 
      height=7)
  print(ggdendrogram(hc_data, rotate = T) +
          theme(axis.text.x = element_text(size = 12, colour = 'black', hjust = .5),
                axis.text.y = element_text(size = 12, colour = 'black'),
                axis.line.x = element_line(size = .5),
                axis.ticks.x = element_line(size = .5),
                axis.title.x = element_text(size = 15, colour = 'black')) +
          labs(y = 'euclidean distance (r)'))
  dev.off()
  
  # multi-dimensional scaling
  fit = cmdscale(dist_mat, eig = T, k = 2)
  x = fit$points[,1]
  y = fit$points[,2]
  mds_plot_data = data.frame(x = x, y = y)
  pdf(file = file.path(resdir, sprintf("MDS_%s.pdf", it)),
      bg="transparent",
      width=7, 
      height=7)
  print(ggplot(mds_plot_data, aes(x=x,y=y))+
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.background = element_blank(),
                panel.border = element_rect(colour = "black", size = 2, fill = "transparent"),
                axis.text.x  = element_blank(),
                axis.text.y  = element_blank(),
                axis.ticks.x  = element_blank(),
                axis.ticks.y  = element_blank(),
                axis.title.x = element_blank(), 
                axis.title.y = element_blank(),
                legend.position="none") +
          geom_point(na.rm=F)+
          geom_text_repel(aes(label = rownames(mds_plot_data), size = 5),box.padding = unit(0.25, "lines")) +
          coord_fixed())
  dev.off()
}

# featquery, face v house in ffa and ppa
data_fq <- data.frame(read.csv(file.path(projDir, 'data', 'mvpa_fq.csv'), header = T))
data_fq_fh <- droplevels(data_fq[data_fq$category %in% c('face','house'),])
data_plot <- aggregate(data = data_fq_fh, signalChange ~ region + imageType + category, FUN = mean)
temp <- aggregate(data = data_fq_fh, signalChange ~ region + imageType + category, FUN = std.error)
data_plot$se <- temp$signalChange
colnames(data_plot)[4] <- 'mean'

resdir <- file.path(projDir, 'results', 'MVPA','fq_face_house')
dir.create(resdir, recursive = T, showWarnings = F)
pdf(file=file.path(resdir, 'means.pdf'),
    bg="transparent",
    width=8, 
    height=6)
print(ggplot(data_plot, aes(x=imageType, y=mean, fill = category)) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line.x = element_line(colour = "black", size = .5), 
              axis.line.y = element_line(colour = "black", size = .5), 
              axis.ticks = element_line(colour = "black", size = .5),
              axis.text.x  = element_text(colour = "black", size=13), 
              axis.text.y  = element_text(colour = "black", size=13), 
              axis.title.x = element_blank(),
              axis.title.y = element_text(colour = 'black', size=15, margin = margin(0,10,0,0)),
              legend.text = element_text(colour = 'black', size = 13),
              legend.title = element_blank(),
              strip.text = element_text(colour = 'black', size = 13)) +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), stat = "identity", width=.3, position=position_dodge(.7), size = .5) +
        geom_bar(stat = "identity", position = "dodge", width = .7, colour = "black") +
        facet_grid(.~region, scales="free", space="free_x", shrink=T) +
        #geom_text(aes(label = label, x = region, y = mean+se+.02), position = position_dodge(width = .7), colour = "black", size = 5) +
        geom_hline(yintercept = 0) +
        scale_fill_manual(values = c('white','grey90')) +
        labs(y = 'signal change (%)'))
#        scale_y_continuous(breaks=seq(ybreaks_min[model], ybreaks_max[model], ybreaks_step[model]), expand=c(0,0)))
dev.off()

# vectorize regions mat and save as csv
inv_cors_mat.avdiag = inv_cors_mat
for (r in 1:16){ for (c in 2:17){ if (r+1 <= c) inv_cors_mat.avdiag[r,c] <- NA}}
cors_vect <- c(na.omit(c(as.matrix(inv_cors_mat.avdiag[,-1]))))
write.csv(cors_vect, '/Volumes/ddc/research/projects/p007/region_matrix_exp2.csv')

