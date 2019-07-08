# p007
# calculates within- versus between-category contrasts for different brain regions

# packages ---------------------------------------------------------------
packages <- c('ggplot2', 'reshape2', 'plyr', 'multcomp', 'broom', 'lsr', 'ez', 'plotrix', 'R.matlab',
              'nlme', 'grid', 'ppcor', 'psych', 'ggrepel', 'ggdendro', 'dendextend', 'car', 'Hmisc')
lapply(packages, require, character.only = TRUE)

# configuration -----------------------------------------------------------
projDir <- '/Volumes/ddc/research/projects/p007/exp1/fMRI'
conditions_names = c("untextured", "textured")

# general objects -------------------------------------------------
conditions = c("grey", "noise")
wang_regions <- c('V1d','V1v','V2d','V2v','V3d','V3v','hV4','VO1','VO2','PHC1','PHC2','V3a','V3b','LO1','LO2','hMT')

clusters = c("01","02","03","04","05","06","07","08","09","10")
clusters_super = c("u01","u02","u03","u04","u05","u06","u07","u08","u09","u10","t01","t02","t03","t04","t05","t06","t07","t08","t09","t10")
grey_cols = c(1:10,21:30,41:50,61:70,81:90,101:110,121:130,141:150,161:170,181:190)
noise_cols = c(211:220,231:240,251:260,271:280,291:300,311:320,331:340,351:360,371:380,391:400)
within_vals = c(1,12,23,34,45,56,67,78,89,100) # values in grey_cols and noise_cols corresponding to within
between_vals = c(11,21:22,31:33,41:44,51:55,61:66,71:77,81:88,91:99) # and for between
between_vals_cluster = matrix(ncol = 10, nrow = 9, # and for between for each cluster 
                              data = c(11,21,31,41,51,61,71,81,91,
                                       11,22,32,42,52,62,72,82,92,
                                       21:22,33,43,53,63,73,83,93,
                                       31:33,44,54,64,74,84,94,
                                       41:44,55,65,75,85,95,
                                       51:55,66,76,86,96,
                                       61:66,77,87,97,
                                       71:77,88,98,
                                       81:88,99,
                                       91:99))
super_cols = c(1:400) # super_matrix
super_within_vals = c(grey_cols[within_vals], noise_cols[within_vals])
super_between_vals = setdiff(super_cols, super_within_vals)

stimdata_names = matrix(data = NA, nrow = 10, ncol = 24)
stimdata_names[1,] = c("trainer", "truck", "laptop", "elbow", "pin", "sunglasses", "fish", "peacock", "spiderweb", "gates", "glove", "springroll", "fence", "lamp", "wolf", "sandcastle", "watch", "toytruck", "gearstick", "pipe", "eagle", "headphones", "lettuce", "mushroom")
stimdata_names[2,] = c("football", "pepper", "tyre", "volleyball", "drainplug", "tape", "coffeecuplid", "baseball", "prawn", "salsa", "wirecoil", "wheel", "cremebrulee", "helmet", "tennisball", "pepper", "dartboard", "die", "discoball", "lime", "die", "clock", "wheel", "stringball")
stimdata_names[3,] = c("brush", "skigoggles", "tapemeasure", "tapemeasure", "rope", "buggy", "truck", "granolabar", "car", "lei", "bracelet", "pasta", "ladybug", "cap", "swimgoggles", "pasta", "rhino", "raisins", "toycar", "towel", "printer", "xylophone", "doorhandle", "shoe")
stimdata_names[4,] = c("tongs", "comb", "shoehorn", "slug", "pliers", "clothespeg", "jackadapter", "secateurs", "boxcutter", "sock", "peapod", "thermometer", "envelope", "boxcutter", "feather", "masher", "crab", "gun", "boxcutter", "tvremote", "battery", "brush", "paperclip", "screwdriver")
stimdata_names[5,] = c("peppershaker", "golfclubs", "spraybottle", "fireextinguisher", "sandtimer", "toybear", "waders", "nutandbolt", "nutandbolt", "statuette", "statue", "mirror", "boiler", "parkingmeter", "mask", "lipstick", "fireextinguisher", "candles", "sacktruck", "trophy", "lipstick", "airpump", "podium", "mobilephone")
stimdata_names[6,] = c("spanner", "garliccrusher", "wok", "ladel", "lollipop", "xylophone", "necklace", "ceilingfan", "banana", "pin", "cucumber", "scissors", "car", "heater", "pen", "swordfish", "crocodile", "toothbrush", "zip", "rail", "spiritlevel", "arm", "screw", "pen")
stimdata_names[7,] = c("toytruck", "toyhorse", "spade", "boot", "highheel", "boot", "showerhead", "ruler", "staircase", "boxcutter", "tapemeasure", "brassknuckles", "acorn", "boot", "lawnmower", "dentalfloss", "pasta", "spade", "teddybear", "owl", "exercisebike", "shoe", "lollipop", "shield")
stimdata_names[8,] = c("turnip", "drinksmachine", "teddybear", "scales", "gastank", "perfume", "lock", "bookend", "cactus", "mailbox", "teddybear", "helmet", "fonduebowl", "shirt", "puppettheatre", "coffeemachine", "capersjar", "shower", "sipcup", "barrell", "clock", "wreath", "statue", "ox")
stimdata_names[9,] = c("holepunch", "pen", "spiritlevel", "legobrick", "shoe", "tern", "legobrick", "babydoll", "comb", "corkscrew", "sunglasses", "umbrella", "spanner", "boxcutter", "highheel", "tvremote", "nutandbolt", "spiritlevel", "plaster", "clothespeg", "narwhal", "ladle", "wheelbarrow", "paintbrush")
stimdata_names[10,] = c("mirror", "polarbear", "headphones", "suitcase", "slidetray", "parachute", "collander", "birdsnest", "watch", "skigoggles", "mask", "hockeyglove", "bottlecap", "brooch", "bun", "padlock", "pencilsharpener", "ladybug", "giftbow", "ribbon", "doughnut", "binoculars", "almonds", "bun")

# functions ---------------------------------------------------------------
sig_code <- function(x){
  y <- c()
  y[x > .05] = ''
  y[x <= .05] = '*'
  y[x <= .01] = '**'
  y[x <= .001] = '***'
  return(y)
}

bootstrap_getP <- function(x){
  bootsamples <- matrix(sample(x, size = 10000 * length(x), replace = TRUE), 10000, length(x))
  bootstats <- apply(bootsamples, 1, mean)
  samplesBelowZero <- bootstats[bootstats<0]
  p <- length(samplesBelowZero)/10000
  return(p)
}

plot_matrix = function(x, label="matrix", unit="value", minmax=signif(c(min(x,na.rm=T), max(x,na.rm=T)), digits = 2), resdir) { # plots matrix
  df = data.frame(cluster = clusters, t(x))
  colnames(df)[2:11] = clusters
  plotdf = melt(df, id = "cluster")
  pdf(file = paste0(file.path(resdir, label), ".pdf"),
      bg="transparent",
      width=4, 
      height=3)
  print(ggplot(plotdf, aes(x=cluster,y=ordered(variable, levels = rev(clusters)),fill=value))+
          geom_tile(na.rm=F)+
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.background = element_blank(),
                axis.text.x  = element_text(colour = "black", size=12),
                axis.text.y  = element_text(colour = "black", size=12),
                axis.title.x = element_text(colour = "black", size=12, margin = margin(10,0,0,0)), 
                axis.title.y = element_text(colour = "black", size=12, margin = margin(0,10,0,0)),
                legend.title = element_text(colour = "black", size=12, angle = 270)) +
          labs(x = "cluster N", y = "cluster N") +
          scale_fill_continuous(low="red",high="yellow", na.value = "white", limits = minmax, breaks=c(minmax), oob = scales::squish,
                                guide = guide_colorbar(nbin=100, ticks=F, barwidth = .5, barheight = 11.4, title.position="right", title = unit, title.hjust=.5))+
          coord_fixed())
  dev.off()
}


# models ----------------------------------------------
# GIST
matrix_gist = matrix(c(.85,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                       -.10,.92,NA,NA,NA,NA,NA,NA,NA,NA,
                       -.45,-.17,.89,NA,NA,NA,NA,NA,NA,NA,
                       .21,-.24,.11,.93,NA,NA,NA,NA,NA,NA,
                       -.02,.01,-.48,-.34,.96,NA,NA,NA,NA,NA, 
                       -.28,-.31,.42,.22,-.30,.95,NA,NA,NA,NA,
                       -.03,-.15,-.23,-.29,-.03,-.17,.90,NA,NA,NA,
                       .09,.40,-.56,-.40,.51,-.51,.03,.88,NA,NA,
                       -.27,-.22,.32,-.24,-.35,.10,.27,-.38,.91,NA,
                       -.36,.21,.53,-.03,-.32,.05,-.33,-.20,.10,.82), ncol = 10, nrow = 10, byrow = T)
vector_gistWB = as.numeric(na.omit(as.vector(t(matrix_gist))))
for (r in 1:10){ for (c in 1:10){ if (r == c){ matrix_gist[r,c] <- NA }}} # option to remove diagonal
vector_gist = as.numeric(na.omit(as.vector(t(matrix_gist))))
vector_gist_NA = as.numeric(as.vector(t(matrix_gist)))
modeldir = '/Volumes/ddc/research/projects/p007/exp1/models'
plot_matrix(matrix_gist, label = "gist", resdir = modeldir, unit = "Pearson's r")
write.csv(vector_gist, file.path(modeldir,'gist_vector.csv'))

# spectral and spatial GISTs
GISTfile2 <- readMat('/Volumes/ddc/research/projects/p007/exp1/stimuli/scripts/clustering/SpectralSpatialGISTMats.mat')
matrix_gistSpect = GISTfile2$matrices[1,,]
for (r in 1:10){ for (c in 1:10){ if (r < c){ matrix_gistSpect[r,c] <- NA }}}
plot_matrix(matrix_gistSpect, label = "gist_Spect", resdir = modeldir, unit = "Pearson's r")
matrix_gistSpat = GISTfile2$matrices[2,,]
for (r in 1:10){ for (c in 1:10){ if (r < c){ matrix_gistSpat[r,c] <- NA }}}
plot_matrix(matrix_gistSpat, label = "gist_Spat", resdir = modeldir, unit = "Pearson's r")

# make plot for supplementary figure, all data hardcoded
plotData <- data.frame(condition = c('untextured', 'textured'), gist = rep(c('spatial\nGIST', 'spectral\nGIST'), each = 2), correlation = c(0.30, 0.38,0.28,0.37))
plotData$condition <- factor(plotData$condition, levels = c('untextured', 'textured'))
pdf(file=(file.path(modeldir, "GISTS_v_neural.pdf")), 
    bg="transparent",
    width=3.5, 
    height=3)
print(ggplot(plotData, aes(x=gist, y=correlation, fill = condition)) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line.x = element_line(colour = "black", size = .5), 
              axis.line.y = element_line(colour = "black", size = .5), 
              axis.ticks = element_line(colour = "black", size = .5),
              axis.text.x  = element_text(colour = "black", size=10), 
              axis.text.y  = element_text(colour = "black", size=10), 
              axis.title.x = element_blank(), 
              axis.title.y = element_text(size=13, margin = margin(0,10,0,0)),
              legend.text = element_text(colour = "black", size=10),
              legend.title = element_blank()) +
        scale_fill_manual(values=c("white", "gray")) +
        geom_bar(aes(fill = condition), stat = "identity", position = "dodge", width = .7, colour = "black") +
        labs(y = "correlation (r)") +
        coord_cartesian(ylim=c(0,.42)) +
        scale_y_continuous(breaks=seq(0,.4,.2), expand = c(0,0)))
dev.off()

# REAL-WORLD SIZE
# code used: 1 = <10cm, 2 = 10-1000cm, 3 = 1000cm+
data_rws = matrix(data = NA, nrow = 10, ncol = 24)
data_rws[1,]=c(2,3,2,2,1,2,2,3,2,3,2,1,2,2,3,2,1,2,2,2,2,2,2,1)
data_rws[2,]=c(2,1,2,2,1,1,1,1,1,1,1,1,2,2,1,1,2,1,2,1,1,2,1,1)
data_rws[3,]=c(2,2,1,1,2,3,3,1,3,2,1,1,1,2,2,1,3,1,1,2,2,2,2,2)
data_rws[4,]=c(2,2,2,1,2,1,1,2,2,2,1,2,2,2,2,2,2,2,2,2,1,2,1,2)
data_rws[5,]=c(2,3,2,2,1,2,3,1,1,2,3,3,3,2,2,1,2,2,3,2,1,2,3,1)
data_rws[6,]=c(2,2,2,2,1,2,2,3,2,1,2,2,3,2,1,3,3,2,1,2,2,2,1,1)
data_rws[7,]=c(2,2,2,2,2,2,2,2,3,2,1,2,1,2,3,1,1,2,2,2,3,2,1,3)
data_rws[8,]=c(1,2,2,2,2,1,1,2,2,3,2,2,2,2,3,2,2,3,1,3,2,2,1,3)
data_rws[9,]=c(2,1,2,1,2,2,1,2,2,1,2,2,2,2,2,2,1,2,1,1,3,2,3,2)
data_rws[10,]=c(2,3,2,2,2,3,2,2,1,2,2,2,1,1,1,1,1,1,1,1,1,2,1,1)

# get number of instances of each size for each cluster
sum_rws = matrix(data=NA, nrow = 10, ncol = 3)
for (size in 1:3){
  sum_rws[,size] = rowSums(data_rws==size)
}

# get dot product and generate matrix
matrix_rws = matrix(data = NA, nrow = 10, ncol = 10)
for (clusA in 1:10){
  for(clusB in 1:10){
    if (clusA >= clusB){
      matrix_rws[clusA,clusB] = sum(sum_rws[clusA,]*sum_rws[clusB,])
    }
  }
}
vector_rws = as.vector(t(matrix_rws))
plot_matrix(matrix_rws, label = "rws", resdir = modeldir, unit = "dot product")
write.csv(vector_rws, file.path(modeldir, 'rws_vector.csv'))


# TOOL/NON-TOOL
# code used: 1 = non-tool, 2 = tool
data_tool = matrix(data = NA, nrow = 10, ncol = 24)
data_tool[1,]=rep(1,24)
data_tool[2,]=rep(1,24)
data_tool[3,]=c(2,1,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
data_tool[4,]=c(2,2,2,1,2,1,1,2,2,1,1,1,1,2,1,2,1,2,2,2,1,2,1,2)
data_tool[5,]=c(2,1,2,2,1,1,1,1,1,1,1,1,1,1,1,1,2,1,2,1,1,2,1,2)
data_tool[6,]=c(2,2,2,2,1,1,1,1,1,1,1,2,1,1,2,1,1,2,1,1,2,1,1,2)
data_tool[7,]=c(1,1,2,1,1,1,1,2,1,2,2,2,1,1,1,1,1,2,1,1,1,1,1,1)
data_tool[8,]=rep(1,24)
data_tool[9,]=c(2,2,2,1,1,1,1,1,2,2,1,2,2,2,1,2,1,2,1,1,1,2,1,2)
data_tool[10,]=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,2,1,1)

# get number of instances of each size for each cluster
sum_tool = matrix(data=NA, nrow = 10, ncol = 2)
for (tool in 1:2){
  sum_tool[,tool] = rowSums(data_tool==tool)
}

# get dot product and generate matrix
matrix_tool = matrix(data = NA, nrow = 10, ncol = 10)
for (clusA in 1:10){
  for(clusB in 1:10){
    if (clusA >= clusB){
      matrix_tool[clusA,clusB] = sum(sum_tool[clusA,]*sum_tool[clusB,])
    }
  }
}
vector_tool = as.vector(t(matrix_tool))
plot_matrix(matrix_tool, label = "tool", resdir = modeldir, unit = "dot product")
write.csv(vector_tool, file.path (modeldir, 'tool_vector.csv'))

# repeat but just sum numbers of tool (not non-tool)
sum_tool2 = rowSums(data_tool==2)

# get dot product and generate matrix
matrix_tool2 = matrix(data = NA, nrow = 10, ncol = 10)
for (clusA in 1:10){
  for(clusB in 1:10){
    if (clusA >= clusB){
      matrix_tool2[clusA,clusB] = sum(sum_tool2[clusA]*sum_tool2[clusB])
    }
  }
}
vector_tool2 = as.vector(t(matrix_tool2))
plot_matrix(matrix_tool2, label = "tool2", resdir = modeldir, unit = "dot product")
write.csv(vector_tool2, file.path(modeldir, 'tool2_vector.csv'))


# SEMANTIC SIMILARITY (WORDNET)
matrix_wordnet_path=t(data.matrix(read.csv(file.path(modeldir,'wordnet_path.csv'), header=F)))
for(x in 1:10){for(y in 1:10){ if (x<y){matrix_wordnet_path[x,y]=NA}}} # remove duplicate values
vector_wordnet_path = as.vector(t(matrix_wordnet_path))
plot_matrix(matrix_wordnet_path, minmax = c(0,1), label = "wordnet_path", resdir = modeldir, unit = "path similarity")

matrix_wordnet_lch=t(data.matrix(read.csv(file.path(modeldir,'wordnet_lch.csv'), header=F)))
for(x in 1:10){for(y in 1:10){ if (x<y){matrix_wordnet_lch[x,y]=NA}}} # remove duplicate values
vector_wordnet_lch = as.vector(t(matrix_wordnet_lch))
plot_matrix(matrix_wordnet_lch, label = "wordnet_lch", resdir = modeldir, unit = "lch similarity")

matrix_wordnet_wup=t(data.matrix(read.csv(file.path(modeldir,'wordnet_wup.csv'), header=F)))
for(x in 1:10){for(y in 1:10){ if (x<y){matrix_wordnet_wup[x,y]=NA}}} # remove duplicate values
vector_wordnet_wup = as.vector(t(matrix_wordnet_wup))
plot_matrix(matrix_wordnet_wup, label = "wordnet_wup", resdir = modeldir, unit = "wup similarity")

# PERCEPTUAL SIMILARITY (CARD-SORTING TASK)
Nsubs=20
data_card = array(data = NA, dim=c(Nsubs,10,10))
for (subj in 1:Nsubs){
  Ldata = data.frame(read.csv(paste0(modeldir,"/behavioural_card/", subj, ".csv"), header = T))
  #  print(levels(Ldata$Label)) # do this to extract pile names for wordle analysis
  for (clusterA in 1:10){
    for (clusterB in 1:10){
      data_card[subj,clusterA,clusterB] = sum(Ldata[,clusterA+3]*Ldata[,clusterB+3])
    }
  }
}
card_ages <- c(62,52,56,56,22,53,41,58,21,20,52,27,20,21,20,20,20,21)

matrix_card = t(apply(data_card, c(2,3), mean, na.rm=T))
for(x in 1:10){for(y in 1:10){ if (x<y){matrix_card[x,y]=NA}}} # remove duplicate values
vector_card = as.vector(t(matrix_card))
plot_matrix(matrix_card, label = "card", resdir = modeldir, unit = "dot product", minmax = c(2.5,14.4))
write.csv(vector_card, file.path(modeldir,'behavioural_card/gmean_vector.csv'))

# MVPA -----------------------------------------------------
for (dataType in c('within')){
  if (dataType == 'within'){
    data <- data.frame(read.table(file.path(projDir, 'data', 'matrices_within.csv'), header = F))
    data <- aggregate(x = data, by = list(data$V1,data$V2,data$V3), FUN = mean)
    data <- data[,c(1:3,8:407)]
  }
  if (dataType == 'between'){
    data <- data.frame(read.table(file.path(projDir, 'data', 'matrices.csv'), header = F))
  }
  colnames(data)[1:3] <- c('subject', 'maskProj', 'region')
  maskProjs = levels(data$maskProj)
  
  for (maskP in maskProjs){
    dataMaskProj=droplevels(data[data$maskProj == maskP,])
    regions = levels(dataMaskProj$region)
    for (region in regions){
      if (region == 'SPL1') next
      dataRegion = dataMaskProj[dataMaskProj$region == region,-(2:3)]
      resdir = file.path(projDir, 'results','MVPA', dataType, maskP, region)
      dir.create(resdir, recursive = TRUE, showWarnings = FALSE)
      
      # neural matrices ----------
      
      # plot correlation matrices for each condition
      mins_maxs = array(NA,dim=c(2,2)) # both matrices need to have same colour scale, get mins and maxs for both
      mean_matrix = colMeans(dataRegion[,2:401], na.rm=T) # get group mean matrices
      mean_matrix[mean_matrix == 'Nan'] <- NA
      for (cond in 1:2){
        cols = get(paste(conditions[cond], "_cols", sep=""))
        temp_mat = t(matrix(mean_matrix[cols], nrow=10, ncol=10))
        assign(paste("matrix",conditions[cond], sep="_"), temp_mat)
        assign(paste("vector", conditions[cond], sep="_"), as.vector(t(temp_mat)))
        write.csv(as.vector(t(temp_mat)), sprintf('%s/%s_vector_r.csv', resdir, conditions[cond]))
        mins_maxs[cond,1]=min(temp_mat, na.rm=T) # maximum and minimum values for colour scales.
        mins_maxs[cond,2]=max(temp_mat, na.rm=T)
      }
      mins_maxs = signif(mins_maxs, digits = 2)
      min_max = signif(c(min(mins_maxs, na.rm=T), max(mins_maxs, na.rm=T)), digits = 2)
      plot_matrix(matrix_grey, minmax=min_max, resdir = resdir, label = "matrix_grey", unit = "Pearson's r")
      plot_matrix(matrix_noise, minmax=min_max, resdir = resdir, label = "matrix_noise", unit = "Pearson's r")
      
      # plot super matrix
      matrix_super = matrix(mean_matrix, ncol=20, nrow=20, byrow = T)
      vector_super = as.vector(t(matrix_super))
      minmax=signif(c(min(matrix_super, na.rm=T),max(matrix_super, na.rm=T)), digits = 2)
      df = data.frame(cluster = clusters_super, t(matrix_super))
      colnames(df)[2:21] = clusters_super
      plotdf = melt(df, id = "cluster")
      plotdf$cluster = factor(plotdf$cluster, levels = clusters_super)
      pdf(file = file.path(resdir, "matrix_super.pdf"),
          bg="transparent",
          width=7, 
          height=6)
      print(ggplot(plotdf, aes(x=cluster,y=ordered(variable, levels = rev(clusters_super)),fill=value))+
              geom_tile(na.rm=F)+
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    axis.text.x  = element_text(colour = "black", size=12, vjust=.5),
                    axis.text.y  = element_text(colour = "black", size=12),
                    axis.title.x = element_text(colour = "black", size=15, margin = margin(10,0,0,0)), 
                    axis.title.y = element_text(colour = "black", size=15, margin = margin(0,10,0,0)),
                    legend.title = element_text(colour = "black", size=15, angle = 270)) +
              labs(x = "untextured                            textured\ncluster N", y = "cluster N\ntextured                            untextured") +
              scale_fill_continuous(low="red",high="yellow", na.value = "white", limits = c(minmax), breaks=c(minmax), oob = scales::squish, guide = guide_colorbar(nbin=100, ticks=F, barwidth = 1, barheight = 25.4, title.position="right", title = "Pearson's r", title.hjust=.5))+
              scale_x_discrete(labels = rep(clusters,2)) +
              scale_y_discrete(labels = rev(rep(clusters,2))) + 
              coord_fixed())
      dev.off()
      
      # pattern reliability ---------------
      
      ## ANALYSES INCLUDING BOTH GREY AND NOISE
      
      # arrange data for stats and bar graphs showing within vs between  
      dataRegionWB <- data.frame(subject = factor(1:length(levels(data$subject))),
                                 grey_within = rowMeans(subset(dataRegion, select = grey_cols[within_vals]+1), na.rm = TRUE),
                                 grey_between = rowMeans(subset(dataRegion, select = grey_cols[between_vals]+1), na.rm = TRUE),
                                 noise_within = rowMeans(subset(dataRegion, select = noise_cols[within_vals]+1), na.rm = TRUE),
                                 noise_between = rowMeans(subset(dataRegion, select = noise_cols[between_vals]+1), na.rm = TRUE))
      dataRegionWB$grey <- dataRegionWB$grey_within-dataRegionWB$grey_between
      dataRegionWB$noise <- dataRegionWB$noise_within-dataRegionWB$noise_between
      dataRegionWB$overall_within <- (dataRegionWB$grey_within+dataRegionWB$noise_within)/2
      dataRegionWB$overall_between <- (dataRegionWB$grey_between+dataRegionWB$noise_between)/2
      dataRegionWB$overall <- dataRegionWB$overall_within-dataRegionWB$overall_between
      
      # next part preps the data for making a bar plot for overall pattern reliability in each condition
      # make template dataframe then add in means and SEs from data
      summary = data.frame(condition = rep(conditions_names, each = 2),
                           comparison = rep(c("within", "between"), 2),
                           mean = c(rep(NA, 4)),
                           se = c(rep(NA, 4)))
      
      summary$mean[1] = mean(dataRegionWB$grey_within); summary$se[1] = std.error(dataRegionWB$grey_within)
      summary$mean[2] = mean(dataRegionWB$grey_between); summary$se[2] = std.error(dataRegionWB$grey_between)
      summary$mean[3] = mean(dataRegionWB$noise_within); summary$se[3] = std.error(dataRegionWB$noise_within)
      summary$mean[4] = mean(dataRegionWB$noise_between); summary$se[4] = std.error(dataRegionWB$noise_between)
      summary$comparison = factor(summary$comparison, levels = c("within", "between"))
      summary$condition = factor(summary$condition, levels = conditions_names)
      
      # create label dataframes to plot significance
      label_data <- data.frame(condition = conditions_names,
                               mean = c((summary$mean[1] + summary$se[1] +.03), 
                                        (summary$mean[3] + summary$se[3] +.03)))
      labels <- c("","","") # significance codes (grey, noise, grey v noise)
      greyWB = t.test(dataRegionWB$grey_within, dataRegionWB$grey_between, paired=T)
      greyWBd = cohensD(dataRegionWB$grey_within, dataRegionWB$grey_between, method='paired')
      if (greyWB$p.value<.05){ labels[1]="*" }
      if (greyWB$p.value<.01){ labels[1]="**" }
      if (greyWB$p.value<.001){ labels[1]="***" }
      noiseWB = t.test(dataRegionWB$noise_within, dataRegionWB$noise_between, paired=T)
      noiseWBd = cohensD(dataRegionWB$noise_within, dataRegionWB$noise_between, method='paired')
      if (noiseWB$p.value<.05){ labels[2]="*" }
      if (noiseWB$p.value<.01){ labels[2]="**" }
      if (noiseWB$p.value<.001){ labels[2]="***" }
      greyNoiseWB = t.test(dataRegionWB$grey, dataRegionWB$noise, paired=T)
      greyNoiseWBd = cohensD(dataRegionWB$grey, dataRegionWB$noise, method='paired')
      if (greyNoiseWB$p.value<.05){ labels[3]="*" }
      if (greyNoiseWB$p.value<.01){ labels[3]="**" }
      if (greyNoiseWB$p.value<.001){ labels[3]="***" }
      label_sizes <- c(5,5,5)
      bar_heights <-c(label_data$mean[1]-.01, label_data$mean[2]-.01, max(label_data$mean)+.03)
      
      # make standard y scale/breaks for HOCSA and Wang
      ylims <- c(-.1, .5)
      ybreaks <- seq(-.1,.5,.1)
      
      # make special limits for P1200 plots
      if (maskP == 'P1200'){
        ylims <- c(-.03,0.21)
        ybreaks <- seq(0,.2,.1)
        label_data <- data.frame(condition = conditions_names,
                                 mean = c((summary$mean[1] + summary$se[1] +.015), 
                                          (summary$mean[3] + summary$se[3] +.015)))
        bar_heights <-c(label_data$mean[1]-.005, label_data$mean[2]-.005, max(label_data$mean)+.015)
      }
      
      # now make plot
      pdf(file=(file.path(resdir, "patt_rel_overall.pdf")), 
          bg="transparent",
          width=3.5, 
          height=3)
      print(ggplot(summary, aes(x=condition, y=mean)) + 
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(), 
                    axis.line.x = element_line(colour = "black", size = .5), 
                    axis.line.y = element_line(colour = "black", size = .5), 
                    axis.ticks = element_line(colour = "black", size = .5),
                    axis.text.x  = element_text(colour = "black", size=10), 
                    axis.text.y  = element_text(colour = "black", size=10), 
                    axis.title.x = element_blank(), 
                    axis.title.y = element_text(size=13, margin = margin(0,10,0,0)),
                    legend.text = element_text(colour = "black", size=10),
                    legend.title = element_blank()) +
              scale_fill_manual(values=c("white", "gray")) +
              geom_errorbar(aes(fill = comparison, ymin=mean-se, ymax=mean+se), stat = "identity", width=.3, position=position_dodge(.7), size = .5) +
              geom_bar(aes(fill = comparison), stat = "identity", position = "dodge", width = .7, colour = "black") +
              geom_text(data = label_data, label = labels[1:2], colour = "black", size = label_sizes[1:2]) +
              geom_text(label = labels[3], x = 1.5, y = bar_heights[3]+.005, colour = "black", size = label_sizes[3]) +
              geom_segment(aes(x=.78, xend=1.22, y=bar_heights[1], yend=bar_heights[1]), size=.5, colour = "black") +
              geom_segment(aes(x=1.78, xend=2.22, y=bar_heights[2], yend=bar_heights[2]), size=.5, colour = "black") +
              geom_segment(aes(x=1, xend=2, y=bar_heights[3], yend=bar_heights[3]), size=.5, colour = "black") +
              geom_hline(yintercept = 0) +
              labs(y = "correlation (z)") +
              coord_cartesian(ylim=ylims) +
              scale_y_continuous(breaks=ybreaks, expand = c(0,0)))
      dev.off()
      
      # make plot for collapse across grey and noise
      summary <- data.frame(comparison = c('within','between'),
                            mean = c(mean(dataRegionWB$overall_within), mean(dataRegionWB$overall_between)),
                            se = c(std.error(dataRegionWB$overall_within), std.error(dataRegionWB$overall_between)))
      summary$comparison <- factor(summary$comparison, levels = c('within','between'))
      
      label <- ''
      overallWB = t.test(dataRegionWB$overall_within, dataRegionWB$overall_between, paired=T)
      overallWBd = cohensD(dataRegionWB$overall_within, dataRegionWB$overall_between, method='paired')
      if (overallWB$p.value<.05){ label="*" }
      if (overallWB$p.value<.01){ label="**" }
      if (overallWB$p.value<.001){ label="***" }
      bar_height <- summary$mean[1]+summary$se[1]+.004
      
      # make standard y scale/breaks for HOCSA and Wang
      ylims <- c(-.1, .5)
      ybreaks <- seq(-.1,.5,.1)
      
      # make special limits for P1200 plots
      if (maskP == 'P1200'){
        ylims <- c(-.05,0.25)
        ybreaks <- seq(0,.2,.1)
      }
      
      # now make plot
      pdf(file=(file.path(resdir, "patt_rel_overall_collapse.pdf")), 
          bg="transparent",
          width=2, 
          height=3)
      print(ggplot(summary, aes(x=comparison, y=mean, fill = rownames(summary))) + 
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(), 
                    axis.line.x = element_line(colour = "black", size = .5), 
                    axis.line.y = element_line(colour = "black", size = .5), 
                    axis.ticks = element_line(colour = "black", size = .5),
                    axis.text.x  = element_text(colour = "black", size=10), 
                    axis.text.y  = element_text(colour = "black", size=10), 
                    axis.title.x = element_blank(), 
                    axis.title.y = element_text(size=13, margin = margin(0,10,0,0)),
                    legend.text = element_text(colour = "black", size=10),
                    legend.title = element_blank()) +
              scale_fill_manual(values=c("white", "gray")) +
              geom_errorbar(aes(ymin=mean-se, ymax=mean+se), stat = "identity", width=.3, position=position_dodge(.7), size = .5) +
              geom_bar(stat = "identity", position = "dodge", width = .7, colour = "black") +
              guides(fill=F) +
              geom_text(label = label, x = 1.5, y = bar_height+.002, colour = "black", size = 5) +
              geom_segment(aes(x=1, xend=2, y=bar_height, yend=bar_height), size=.5, colour = "black") +
              geom_hline(yintercept = 0) +
              labs(y = "correlation (z)") +
              coord_cartesian(ylim=ylims) +
              scale_y_continuous(breaks=ybreaks, expand=c(0,0)))
      dev.off()
      
      # arrange data for ANOVA, with condition, cluster and comparison as repeated measures
      # calculate within and between for each cluster
      for (cluster in 1:10){
        dataRegionWB[,(10+(cluster*4)-3)] = dataRegion[,grey_cols[within_vals[cluster]]+1]
        dataRegionWB[,(11+(cluster*4)-3)] = rowMeans(dataRegion[,grey_cols[between_vals_cluster[,cluster]]+1])
        dataRegionWB[,(12+(cluster*4)-3)] = dataRegion[,noise_cols[within_vals[cluster]]+1]
        dataRegionWB[,(13+(cluster*4)-3)] = rowMeans(dataRegion[,noise_cols[between_vals_cluster[,cluster]]+1])
        colnames(dataRegionWB)[(10+(cluster*4)-3)] = paste0("grey_cl", cluster, "_w")
        colnames(dataRegionWB)[(11+(cluster*4)-3)] = paste0("grey_cl", cluster, "_b")
        colnames(dataRegionWB)[(12+(cluster*4)-3)] = paste0("noise_cl", cluster, "_w")
        colnames(dataRegionWB)[(13+(cluster*4)-3)] = paste0("noise_cl", cluster, "_b")
      }
      
      # lengthen data frame, z-transform and run ANOVA
      temp=dataRegionWB[,-(2:10)]
      dataRegionWB_long = melt(temp, id.vars = "subject")
      colnames(dataRegionWB_long)[2:3] = c("condition", "r")
      dataRegionWB_long$z = fisherz(dataRegionWB_long$r)
      dataRegionWB_long$bg = factor(rep(conditions, each = length(levels(data$subject))*2))
      dataRegionWB_long$comparison = factor(rep(c("within", "between"), each = length(levels(data$subject))))
      dataRegionWB_long$cluster = factor(rep(clusters, each = length(levels(data$subject))*4))
      dataRegionWB_long$condID = factor(with(dataRegionWB_long, paste0(comparison, cluster)))
      ezmodel <- ezANOVA(data = dataRegionWB_long, dv = z, wid = subject, within = .(bg, comparison, cluster), detailed = TRUE)
      
      
      # write results to text file
      anovaFile <- file.path(resdir, "ANOVA.txt")
      sink(anovaFile)
      cat(sprintf('### ANOVA ###\n'))
      cat(capture.output(ezmodel), sep = '\n')
      cat(sprintf('# overall W>B for grey and noise #\n'))
      cat(capture.output(overallWB), sep = '\n')
      cat(sprintf('Cohens D: %s\n', capture.output(overallWBd)), sep = '\n')
      cat(sprintf('# grey W>B #\n'))
      cat(capture.output(greyWB), sep = '\n')
      cat(sprintf('Cohens D: %s\n', capture.output(greyWBd)), sep = '\n')
      cat(sprintf('# noise W>B #\n'))
      cat(capture.output(noiseWB), sep = '\n')
      cat(sprintf('Cohens D: %s\n', capture.output(noiseWBd)), sep = '\n')
      cat(sprintf('# grey vs noise W>B #\n'))
      cat(capture.output(greyNoiseWB), sep = '\n')
      cat(sprintf('Cohens D: %s\n', capture.output(greyNoiseWBd)), sep = '\n')
      sink()
      
      ## SEPARATE ANALYSES FOR GREY AND NOISE
      
      # set up pw_comps to see which clusters are significantly distinct
      for (cond in conditions){ # for each scan condition
        comp_df <- data.frame(cluster = 1:10, mean = NA, se = NA, t = NA, df = 20, CI_lo = NA, CI_hi = NA, p = NA, p.sig = NA, p.cor = NA, p.cor.sig = NA, cohensD = NA)
        dataCond <- droplevels(dataRegionWB_long[dataRegionWB_long$bg == cond,]) # restrict data to this condition
        pw <- pairwise.t.test(dataCond$z, dataCond$condID, p.adjust.method = 'none')
        for (c in 1:10){
          within <- dataCond$z[dataCond$cluster == sprintf("%02d", c) & dataCond$comparison == 'within']
          between <- dataCond$z[dataCond$cluster == sprintf("%02d", c) & dataCond$comparison == 'between']
          test <- t.test(within, between, paired = T)
          wb <- within - between
          comp_df[c,-1] <- c(mean(wb), std.error(wb), test$statistic, 20, test$conf.int[1], test$conf.int[2], pw$p.value[sprintf('within%02d', c), sprintf('between%02d', c)], NA, NA, NA, cohensD(within, between, method = 'pooled'))
        }
        comp_df$p.sig <- sig_code(comp_df$p) # add significance codes
        comp_df$p.cor <- p.adjust(comp_df$p, method = 'holm') # calculated corrected p values
        comp_df$p.cor.sig  <- sig_code(comp_df$p.cor)
        
        # add results to ANOVA text file
        sink(anovaFile, append = T)
        cat(sprintf('\n# pairwise contrast - %s #\n', cond))
        cat(capture.output(comp_df), sep = '\n')
        sink() # unsink text file
        
        # now do permutation testing version
        temp_cols <- get(sprintf('%s_cols', cond))
        permutations <- array(data = NA, dim = c(length(levels(data$subject)),10,1000))
        originals <- matrix(data = NA, nrow = 10, ncol = length(levels(data$subject)))
        withinsAD <- c(1,3,6,10,15,21,28,36,45,55)
        betweensAD <- matrix(data = c(2,4,7,11,16,22,29,37,46,
                                      2,5,8,12,17,23,30,38,47,
                                      4,5,9,13,18,24,31,39,48,
                                      7:9,14,19,25,32,40,49,
                                      11:14,20,26,33,41,50,
                                      16:20,27,34,42,51,
                                      22:27,35,43,52,
                                      29:35,44,53,
                                      37:44,54,
                                      46:54), nrow = 10, ncol = 9, byrow = T)
        
        for (s in 1:length(levels(data$subject))){
          temp_mat <- fisherz(as.numeric(dataRegion[s,temp_cols+1]))
          temp_mat <- temp_mat[!is.na(temp_mat)]
          for (c in 1:10){
            originals[c,s] <- temp_mat[withinsAD[c]] - mean(temp_mat[betweensAD[c,]])
            for (p in 1:1000){
              shuf_mat <- sample(temp_mat, size = length(temp_mat), replace = F)
              permutations[s,c,p] <- shuf_mat[withinsAD[c]] - mean(shuf_mat[betweensAD[c,]])
            }
          }
        }
        permutations.collapse <- apply(permutations, c(2,3), mean)
        originals.collapse <- apply(originals, 1, mean)
        dataPerm <- data.frame(cluster = 1:10, p = NA, p.sig = " ", p.cor = NA, p.cor.sig = " ")
        for (c in 1:10){
          dataPerm$p[c] <- (length(permutations.collapse[c,permutations.collapse[c,]>originals.collapse[c]])+1)/1000 # add 1 to include original matrix
        }
        dataPerm$p.sig <- sig_code(dataPerm$p)
        dataPerm$p.cor <- p.adjust(dataPerm$p, method = 'fdr')
        dataPerm$p.cor.sig <- sig_code(dataPerm$p.cor)
        
        # make summary data frame for plotting
        plot_data = data.frame(cluster = clusters, mean = comp_df$mean, se = comp_df$se, sig = comp_df$p.cor.sig)
        pdf(file=sprintf('%s/patt_rel_cluster_%s.pdf', resdir, cond), 
            bg="transparent",
            width=4, 
            height=3)
        print(ggplot(plot_data, aes(x=cluster, y=mean)) + 
                theme(panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(), 
                      panel.background = element_blank(), 
                      axis.line.x = element_line(colour = "black", size = .5), 
                      axis.line.y = element_line(colour = "black", size = .5), 
                      axis.ticks = element_line(colour = "black", size = .5),
                      axis.text.x  = element_text(colour = "black", size=10), 
                      axis.text.y  = element_text(colour = "black", size=10), 
                      axis.title.x = element_text(size=13, margin = margin(10,0,0,0)), 
                      axis.title.y = element_text(size=13, margin = margin(0,10,0,0))) +
                geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.3, position=position_dodge(.9), size = .5) +
                geom_bar(stat = "identity", position = "dodge", width = .7, fill = "gray90", colour = "black") +
                geom_text(aes(label = sig, y = mean+se+.03), colour = "black", size = 5) +
                geom_hline(yintercept = 0) +
                labs(x = "cluster N", y = "within minus between (z)") +
                scale_y_continuous(expand = c(0, 0), limits = c(-.1,.5), breaks=c(0,.2,.4,.6,.8))) # limits = c(0,.85),
        dev.off()
        
        
        # add results to ANOVA text file
        sink(anovaFile, append = T)
        cat(sprintf('\n# permutation test for each cluster - %s #\n', cond))
        cat(capture.output(dataPerm), sep = '\n')
        sink() # unsink text file
      }
      
      # intermatrix correlations
      nmats = 12
      # put neural and model matrices into matrix
      matrices = matrix(data = NA, nrow = nmats, ncol = 100)
      matrices[1,] <- fisherz(vector_grey) # z transform matrices containing R values
      matrices[2,] <- fisherz(vector_noise)
      matrices[3,] <- fisherz(vector_gist_NA)
      matrices[4,] <- vector_card
      matrices[5,] <- vector_wordnet_path
      matrices[6,] <- vector_wordnet_lch
      matrices[7,] <- vector_wordnet_wup
      matrices[8,] <- vector_rws
      matrices[9,] <- vector_tool
      matrices[10,] <- vector_tool2
      matrices[11,] <- c(t(matrix_gistSpat))
      matrices[12,] <- c(t(matrix_gistSpect))
      
      # save out between-value only matrices
      write.csv(matrices[1,between_vals], sprintf('%s/grey_vector_between_z.csv', resdir))
      write.csv(matrices[2,between_vals], sprintf('%s/noise_vector_between_z.csv', resdir))
      
      # mean-centre the within and between separately to control for W>B in matrix correlations
      matrices_demean <- matrices
      matrices_demean[,within_vals] <- matrices_demean[,within_vals]-rowMeans(matrices_demean[,within_vals])
      matrices_demean[,between_vals] <- matrices_demean[,between_vals]-rowMeans(matrices_demean[,between_vals])
      
      # save out demeaned matrices
      write.csv(matrices_demean[1,], sprintf('%s/grey_vector_demean_z.csv', resdir))
      write.csv(matrices_demean[2,], sprintf('%s/noise_vector_demean_z.csv', resdir))
      
      # remove NAs
      matrices_noNA <- matrix(matrices[!is.na(matrices)],nmats,55)
      matrices_demean_noNA <- matrix(matrices_demean[!is.na(matrices_demean)],nmats,55)
      
      # get correlations between matrices
      library(Hmisc)
      mat_cor_tests <- Hmisc::rcorr(t(matrices_noNA))
      mat_cor_tests_demean <- Hmisc::rcorr(t(matrices_demean_noNA))
      
      # label the matrix and save out in text file
      matrix_names = c("untextured", "textured", "GIST", "perceptual", "semantic", "wordnet_lch", "wordnet_wup", "rws", "tool", "tool2", 'gist_spatial','gist_spectral')
      mat_cor_coefs <- data.frame(matrixID = matrix_names, mat_cor_tests$r)
      mat_cor_pvals <- data.frame(matrixID = matrix_names, mat_cor_tests$P)
      colnames(mat_cor_coefs) <- c("matrixID", matrix_names)
      colnames(mat_cor_pvals) <- c("matrixID", matrix_names)
      sink(file.path(resdir, "matrix_cors.txt"))
      cat(sprintf('### matrix correlations (within and between) ###\n'))
      cat(capture.output(mat_cor_coefs), sep = '\n')
      cat(sprintf('\n# matrix p values #\n'))
      cat(capture.output(mat_cor_pvals), sep = '\n')
      sink()
      
      # repeat but leave out within-category correlations
      mat_cor_tests <- rcorr(t(matrices[,between_vals]))
      mat_cor_coefs <- data.frame(matrixID = matrix_names, mat_cor_tests$r)
      mat_cor_pvals <- data.frame(matrixID = matrix_names, mat_cor_tests$P)
      colnames(mat_cor_coefs) <- c("matrixID", matrix_names)
      colnames(mat_cor_pvals) <- c("matrixID", matrix_names)
      sink(file.path(resdir, "matrix_cors.txt"), append = T)
      cat(sprintf('\n\n### matrix correlations (between only) ###\n'))
      cat(capture.output(mat_cor_coefs), sep = '\n')
      cat(sprintf('\n# matrix p values #\n'))
      cat(capture.output(mat_cor_pvals), sep = '\n')
      sink()
      
      # label the demeaned matrix and save out in text file
      mat_cor_coefs <- data.frame(matrixID = matrix_names, mat_cor_tests_demean$r)
      mat_cor_pvals <- data.frame(matrixID = matrix_names, mat_cor_tests_demean$P)
      colnames(mat_cor_coefs) <- c("matrixID", matrix_names)
      colnames(mat_cor_pvals) <- c("matrixID", matrix_names)
      sink(file.path(resdir, "matrix_cors.txt"), append = T)
      cat(sprintf('\n\n### matrix correlations (within and between, demeaned) ###\n'))
      cat(capture.output(mat_cor_coefs), sep = '\n')
      cat(sprintf('\n# matrix p values #\n'))
      cat(capture.output(mat_cor_pvals), sep = '\n')
      sink()
      
      # plot matrix correlations for key matrices
      matrices2use = (1:5)
      plot_df = mat_cor_coefs[matrices2use,c(1,matrices2use+1)]
      # remove duplicate values (make into triangle)
      for (r in 1:length(matrices2use)){
        for (c in 1:length(matrices2use)){
          if (r > c){
            plot_df[r,c+1] = NA
          }
        }
      }
      plot_data = melt(plot_df, id = "matrixID")
      plot_data$matrixID = factor(plot_data$matrixID, levels = matrix_names)
      plot_data$variable = factor(plot_data$variable, levels = rev(matrix_names))
      plot_data$label <-  ''
      plot_data$label[plot_data$value > 0.2656] <- '*' # p < .05
      plot_data$label[plot_data$value > 0.3446] <- '**' # p < .01
      plot_data$label[plot_data$value > 0.4317] <- '***' # p < .001
      plot_data$label[plot_data$value == 1] <- '' # remove text for diagonal
      minmax = round(c(min(plot_data$value, na.rm = T), max(plot_data$value, na.rm = T)), digits = 2)
      pdf(file = file.path(resdir, 'matrix_cors.pdf'), bg="transparent", width=4, height=3)
      print(ggplot(plot_data, aes(x=matrixID,y=variable,fill=value)) +
              geom_tile(na.rm=F)+
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(), 
                    panel.background = element_blank(),
                    axis.text.x  = element_text(colour = "black", size=12, angle = 45, hjust = 1, vjust = 1),
                    axis.text.y  = element_text(colour = "black", size=12),
                    axis.title.x = element_blank(), 
                    axis.title.y = element_blank(),
                    legend.title = element_text(colour = "black", size=12, angle = 270)) +
              scale_fill_continuous(low="red",high="yellow", na.value = "white", limits = minmax, breaks=minmax, oob = scales::squish,
                                    guide = guide_colorbar(nbin=100, ticks=F, barwidth = .5, barheight = 9.7, title.position="right", title = "Pearson's r", title.hjust=.5)) +
              geom_text(aes(label = label), nudge_y = -.1) + # significance labels
              geom_hline(yintercept = 0) +
              coord_fixed())
      dev.off()
      
      # run partial correlations
      matrices2use <- (1:5)
      partialCorData = as.data.frame(t(matrices_noNA[matrices2use,]))
      colnames(partialCorData) <- matrix_names[matrices2use]
      library(gtools)
      contrasts <- permutations(n=5,r=3,v=(1:5), repeats.allowed=F)
      
      sink(file.path(resdir, "matrix_cors.txt"), append = T)
      cat('\n### partial correlations ###\n')
      sink()
      
      for (contrast in 1:dim(contrasts)[1]){
        tempMats <- contrasts[contrast,]
        # run partial correlation
        partialCor = pcor.test(partialCorData[,tempMats[1]], partialCorData[,tempMats[2]], partialCorData[,tempMats[3]])
        # write to correlations txt file
        sink(file.path(resdir, "matrix_cors.txt"), append = T)
        cat(sprintf('\n# Correlation between %s and %s, controlling for %s #\n', matrix_names[tempMats[1]], matrix_names[tempMats[2]], matrix_names[tempMats[3]]))
        cat(capture.output(partialCor), sep = '\n')
        sink()
      }
      
      # make scatterplots between untextured, textured, GIST, card and semantic
      titles = c("fMRI response - untextured (z)", "fMRI response - textured (z)", "GIST similarity (z)", "perceptual similarity (dot product)", "semantic similarity (arbitrary units)")
      for (cond in 1:2){
        for (model in 3:5){
          plot_df = data.frame(matrixA = matrices[cond,between_vals], matrixB = matrices[model,between_vals])
          contrast_name = paste(matrix_names[cond], matrix_names[model], sep = "_")
          pdf(file=paste0(resdir, "/scatter_", contrast_name, ".pdf"), 
              bg="transparent",
              width=5, 
              height=5)
          print(ggplot(plot_df, aes(x=matrixA, y=matrixB)) + 
                  geom_smooth(method = "lm", alpha = .15, size = 1, fill = 'dodgerblue3', col='dodgerblue3') + 
                  geom_point(size = 3) +
                  labs(x = conditions[cond], y = matrix_names[model]) + 
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
                  labs(x = titles[cond], y = titles[model]))
          dev.off()
        }
      }
    }
  }
}


# make matrix showing correlations between neural and model matrices
matrix_names_RegsMods <- c('VVP', 'FFA', 'PPA', 'GIST', 'perceptual')
matrices_RegsMods <- matrix(data = c(read.csv(file.path(projDir, 'results/MVPA/between/HOCSA/ventral_stream/grey_vector_between_z.csv'))$x,
                                 read.csv(file.path(projDir, 'results/MVPA/between/P1200/ffa_bi_300/grey_vector_between_z.csv'))$x,
                                 read.csv(file.path(projDir, 'results/MVPA/between/P1200/ppa_bi_300/grey_vector_between_z.csv'))$x,
                                 matrices[3,between_vals],
                                 matrices[4,between_vals]), nrow = 45, ncol = 5)
cor_tests_RegsMods <- rcorr(matrices_RegsMods, type = 'pearson')
cor_coefs_RegsMods <- data.frame(matrixID = matrix_names_RegsMods, cor_tests_RegsMods$r)
colnames(cor_coefs_RegsMods)[2:6] <- matrix_names_RegsMods
for (r in 1:dim(matrices_RegsMods)[2]){for (c in 1:dim(matrices_RegsMods)[2]){if (r >= c){cor_coefs_RegsMods[r,c+1] = NA}}} # remove duplicate values (make into triangle)
plot_data <- melt(cor_coefs_RegsMods, id = "matrixID")
plot_data$matrixID <- factor(plot_data$matrixID, levels = matrix_names_RegsMods)
plot_data$variable <- factor(plot_data$variable, levels = rev(matrix_names_RegsMods))
plot_data$label <-  ''
plot_data$label[plot_data$value > 0.2656] <- '*' # p < .05
plot_data$label[plot_data$value > 0.3446] <- '**' # p < .01
plot_data$label[plot_data$value > 0.4317] <- '***' # p < .001
plot_data$label[plot_data$value == 1] <- '' # remove text for diagonal
minmax = round(c(min(plot_data$value, na.rm = T), max(plot_data$value, na.rm = T)), digits = 2)
pdf(file = file.path(projDir, 'matrix_fig.pdf'), bg="transparent", width=4, height=3)
print(ggplot(plot_data, aes(x=matrixID,y=variable,fill=value)) +
        geom_tile(na.rm=F)+
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(),
              axis.text.x  = element_text(colour = "black", size=12, angle = 45, hjust = 1, vjust = 1),
              axis.text.y  = element_text(colour = "black", size=12),
              axis.title.x = element_blank(), 
              axis.title.y = element_blank(),
              legend.title = element_text(colour = "black", size=12, angle = 270)) +
        scale_fill_continuous(low="red",high="yellow", na.value = "white", limits = minmax, breaks=minmax, oob = scales::squish,
                              guide = guide_colorbar(nbin=100, ticks=F, barwidth = .5, barheight = 9.7, title.position="right", title = "Pearson's r", title.hjust=.5)) +
        geom_text(aes(label = label), nudge_y = -.1) + # significance labels
        coord_fixed())
dev.off()



# W v B for each model ---------------------------------------------------

models <- c('GIST', 'perceptual', 'semantic')

# replace the diagonal (removed earlier)
matrix_gist = matrix(c(.85,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                       -.10,.92,NA,NA,NA,NA,NA,NA,NA,NA,
                       -.45,-.17,.89,NA,NA,NA,NA,NA,NA,NA,
                       .21,-.24,.11,.93,NA,NA,NA,NA,NA,NA,
                       -.02,.01,-.48,-.34,.96,NA,NA,NA,NA,NA, 
                       -.28,-.31,.42,.22,-.30,.95,NA,NA,NA,NA,
                       -.03,-.15,-.23,-.29,-.03,-.17,.90,NA,NA,NA,
                       .09,.40,-.56,-.40,.51,-.51,.03,.88,NA,NA,
                       -.27,-.22,.32,-.24,-.35,.10,.27,-.38,.91,NA,
                       -.36,.21,.53,-.03,-.32,.05,-.33,-.20,.10,.82), ncol = 10, nrow = 10, byrow = T)

# put model matrices into 3D matrix
model_mats <- array(NA, c(3,10,10))
model_mats[1,,] <- matrix_gist
model_mats[2,,] <- matrix_card
model_mats[3,,] <- matrix_wordnet_path

# get matrix coordinates for each cluster's between values
between_vals_r = matrix(nrow=10,ncol=9,data=c(2,3,4,5,6,7,8,9,10,
                                              2,3,4,5,6,7,8,9,10,
                                              3,3,4,5,6,7,8,9,10,
                                              4,4,4,5,6,7,8,9,10,
                                              5,5,5,5,6,7,8,9,10,
                                              6,6,6,6,6,7,8,9,10,
                                              7,7,7,7,7,7,8,9,10,
                                              8,8,8,8,8,8,8,9,10,
                                              9,9,9,9,9,9,9,9,10,
                                              10,10,10,10,10,10,10,10,10))
between_vals_c = matrix(nrow=10,ncol=9,data=c(1,1,1,1,1,1,1,1,1,
                                              1,2,2,2,2,2,2,2,2,
                                              1,2,3,3,3,3,3,3,3,
                                              1,2,3,4,4,4,4,4,4,
                                              1,2,3,4,5,5,5,5,5,
                                              1,2,3,4,5,6,6,6,6,
                                              1,2,3,4,5,6,7,7,7,
                                              1,2,3,4,5,6,7,8,8,
                                              1,2,3,4,5,6,7,8,9,
                                              1,2,3,4,5,6,7,8,9))

# set up blank dataframe in which to put within and average between cluster values for each cluster * model
model_WB <- data.frame(model = rep(models, each = 20),
                       cluster = rep(1:10, each = 2),
                       comparison = c('within', 'between'),
                       value = NA)

# input values
for (model in 1:3){
  for (cluster in 1:10){
    model_WB$value[((model-1)*20)+((cluster*2)-1)] <- model_mats[model,cluster,cluster] # within vals
    betweens <- c()
    for (comp in 1:9){
      betweens[comp] <- model_mats[model,between_vals_r[cluster,comp],between_vals_c[cluster,comp]]
      model_WB$value[((model-1)*20)+(cluster*2)] <- mean(betweens) # between vals
    }
  }
}

# run t tests, output results, make plot
sink(file.path(modeldir,'model_WB.txt'))
cat('### t-tests for W v B for each model ###\n')
sink()

titles = c("fMRI response - untextured (z)", "fMRI response - textured (z)", "GIST similarity (z)", "perceptual similarity (dot product)", "semantic similarity (arbitrary units)")
titles_y = titles[3:5]

# limits and breaks of y axis for each model
ylims_min <- c(-.2,0,0)
ylims_max <- c(1.1,13,0.13)
ybreaks_min <- c(0,0,0)
ybreaks_max <- c(1,12,0.12)
ybreaks_step <- c(.5,4,.04)

for (model in 1:3){
  
  # t tests
  model_data <- droplevels(model_WB[model_WB$model == models[model],])
  model_test <- t.test(model_data$value[model_data$comparison == 'within'], model_data$value[model_data$comparison == 'between'], paired = T)
  model_d = cohensD(model_data$value[model_data$comparison == 'within'], model_data$value[model_data$comparison == 'between'], method='paired')
  
  # output results
  sink(file.path(modeldir,'model_WB.txt'), append = T)
  cat(sprintf('\n# %s #\n', models[model]))
  cat(capture.output(model_test), sep = '\n')
  cat(sprintf('Cohens D: %s\n', capture.output(model_d)), sep = '\n')
  sink()
  
  # make plot
  plot_df <- aggregate(value ~ comparison, data = model_data, FUN = mean)
  temp <- aggregate(value ~ comparison, data = model_data, FUN = std.error)
  plot_df$se <- temp$value
  colnames(plot_df)[2] <- 'mean'
  plot_df$comparison <- factor(plot_df$comparison, levels = c('within', 'between'))
  bar_height <- max(plot_df$mean)+plot_df$se[plot_df$mean == max(plot_df$mean)] + (ylims_max[model]-ylims_min[model])/30 # significance bar height is largest value plus error bar plus 1/20th of the data range
  label_height <- bar_height + (ylims_max[model]-ylims_min[model])/40
  label <- ''
  if (model_test$p.value<.05){ label <- "*" }
  if (model_test$p.value<.01){ label <- "**" }
  if (model_test$p.value<.001){ label <- "***" }
  ybreaks_temp <- seq
  
  pdf(file=(sprintf('%s/%s_WB.pdf', modeldir, models[model])), 
      bg="transparent",
      width=2, 
      height=3)
  print(ggplot(plot_df, aes(x=comparison, y=mean, fill = rownames(plot_df))) + 
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.background = element_blank(), 
                axis.line.x = element_line(colour = "black", size = .5), 
                axis.line.y = element_line(colour = "black", size = .5), 
                axis.ticks = element_line(colour = "black", size = .5),
                axis.text.x  = element_text(colour = "black", size=10), 
                axis.text.y  = element_text(colour = "black", size=10), 
                axis.title.x = element_blank(), 
                axis.title.y = element_text(size=13, margin = margin(0,10,0,0)),
                legend.text = element_text(colour = "black", size=10),
                legend.title = element_blank()) +
          scale_fill_manual(values=c("gray", "white")) +
          geom_errorbar(aes(ymin=mean-se, ymax=mean+se), stat = "identity", width=.3, position=position_dodge(.7), size = .5) +
          geom_bar(stat = "identity", position = "dodge", width = .7, colour = "black") +
          guides(fill=F) +
          geom_text(label = label, x = 1.5, y = label_height, colour = "black", size = 5) +
          geom_segment(aes(x=1, xend=2, y=bar_height, yend=bar_height), size=.5, colour = "black") +
          labs(y = titles_y[model]) +
          geom_hline(yintercept = 0) +
          coord_cartesian(ylim=c(ylims_min[model], ylims_max[model])) +
          scale_y_continuous(breaks=seq(ybreaks_min[model], ybreaks_max[model], ybreaks_step[model]), expand=c(0,0)))
  dev.off()
}




# Wang regions - W v B, region matrix and hierarchical clustering --------------------------------------------------

for (cond in 1:2){
  if (cond == 1) { cols2use <- grey_cols; bg <- 'grey' }
  if (cond == 2) { cols2use <- noise_cols; bg <- 'noise' }
  
  # W v B
  WB_wang_data <- data.frame(region = wang_regions,
                           stream = c(rep('posterior',7), rep('ventral',4), rep('lateral',5)),
                           mean = NA, se = NA, p = NA, t = NA, p.sig = NA, p.cor = NA, p.cor.sig = NA, CohensD = NA)
  
  for (r in 1:dim(WB_wang_data)[1]){
    temp_data <- data[data$region == as.character(WB_wang_data$region[r]), cols2use+3]
    within <- fisherz(temp_data[,within_vals])
    within_means <- rowMeans(within, na.rm=T)
    between <- fisherz(temp_data[,between_vals])
    between_means <- rowMeans(between, na.rm=T)
    WB_means <- within_means-between_means
    WB_wang_data$mean[r] <- mean(WB_means)
    WB_wang_data$se[r] <- std.error(WB_means)
    ttest <- t.test(WB_means, mu = 0)
    WB_wang_data$p[r] <- ttest$p.value
    WB_wang_data$t[r] <- ttest$statistic
    WB_wang_data$CohensD[r] <- cohensD(within_means, between_means, method='paired')
  }
  WB_wang_data$p.cor <- p.adjust(c(WB_wang_data$p), method = 'holm', n = length(WB_wang_data$p)) # correct p values
  WB_wang_data$p.sig <- sig_code(WB_wang_data$p)
  WB_wang_data$p.cor.sig <- sig_code(WB_wang_data$p.cor)
  
  # put results in text file
  sink(sprintf('%s/results/MVPA/between/Wang_2015/all_regions/WB_%s.txt', projDir, bg))
  cat('### t-tests for W v B for each region ###\n')
  cat(capture.output(WB_wang_data), sep = '\n')
  sink()
  
  WB_wang_data$region <- factor(WB_wang_data$region, levels = wang_regions)
  WB_wang_data$stream <- factor(WB_wang_data$stream, levels = c('posterior','ventral','lateral'))
  
  pdf(file=sprintf('%s/results/MVPA/between/Wang_2015/all_regions/WB_all_regions_%s.pdf', projDir, bg),
      bg="transparent",
      width=10, 
      height=3)
  print(ggplot(WB_wang_data, aes(x=region, y=mean)) + 
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.background = element_blank(), 
                axis.line.x = element_line(colour = "black", size = .5), 
                axis.line.y = element_line(colour = "black", size = .5), 
                axis.ticks = element_line(colour = "black", size = .5),
                axis.text.x  = element_text(colour = "black", size=10), 
                axis.text.y  = element_text(colour = "black", size=10), 
                axis.title.x = element_blank(),
                axis.title.y = element_text(size=13, margin = margin(0,10,0,0))) +
          geom_errorbar(aes(ymin=mean-se, ymax=mean+se), stat = "identity", width=.3, position=position_dodge(.7), size = .5) +
          geom_bar(stat = "identity", position = "dodge", width = .7, colour = "black", fill = "white") +
          facet_grid(.~stream, scales="free", space="free_x", shrink=T) +
          geom_text(aes(label = label, y = mean+se+.02), colour = "black", size = 5) +
          geom_hline(yintercept = 0) +
          labs(y = 'within minus between (z)'))
  #        scale_y_continuous(breaks=seq(ybreaks_min[model], ybreaks_max[model], ybreaks_step[model]), expand=c(0,0)))
  dev.off()
}

# region matrix
Nregions = length(wang_regions)
all_mean_mats = matrix(data = NA, nrow = Nregions, ncol = 45) # set up blank matrix
rownames(all_mean_mats) = wang_regions
for (reg in 1:Nregions){
  temp_mat <- data.frame(read.csv(sprintf('%s/results/MVPA/between/Wang_2015/%s/grey_vector_between_z.csv', projDir, wang_regions[reg]), header = T))
  all_mean_mats[reg,] <- temp_mat$x
}
inv_cors_mat <- data.frame(region = wang_regions, 1-(cor(t(all_mean_mats)))) # matrix of 1-R
inv_cors_mat[inv_cors_mat == 0] = NA
temp_plot_data <- melt(inv_cors_mat, id = "region")
temp_plot_data$region <- factor(temp_plot_data$region, levels = wang_regions) # change order of regions
temp_plot_data$variable <- factor(temp_plot_data$variable, levels = rev(wang_regions))
minmax=c(min(temp_plot_data$value, na.rm=T)-.01, max(temp_plot_data$value, na.rm=T)+.01)
rounded_minmax <- round(minmax, digits=2)
pdf(file = file.path(projDir, "results/MVPA/between/Wang_2015/all_regions/region_matrix.pdf"),
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
        scale_fill_continuous(low="yellow",high="red", na.value = "grey", limits = c(0.03,1.15), breaks=c(0.03,1.15),
                              guide = guide_colorbar(nbin=100, ticks=F, barwidth = 1, barheight = 26, title.position="right", title = "distance (1-r)", title.hjust=.5))+
        coord_fixed())
dev.off()

# hierarchical clustering
dist_mat <- dist(all_mean_mats, method = 'maximum') # matrix of euclidean distances
hc_data <- hclust(dist_mat, method = 'complete')
pdf(file = file.path(projDir,"results/MVPA/between/Wang_2015/all_regions/hier_clust.pdf"),
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


# card figure -------------------------------------------------------------
subj2use <- 6
Ldata = data.frame(read.csv(paste0(modeldir, "/behavioural_card/", subj2use, ".csv"), header = T))
df = data.frame(cluster = clusters, t(Ldata[,4:13]))
colnames(df)[2:11] <- as.character(Ldata$Label)
plotdf = melt(df, id = "cluster")
plotdf$value = cut(plotdf$value, breaks=6, labels=0:5, include.lowest = T)
plotdf$value <- factor(plotdf$value, levels = rev(levels(plotdf$value)))
pdf(file = file.path(projDir, "report/figures/wordle/matrix_count.pdf"),
    bg="transparent",
    width=5, 
    height=3.5)
print(ggplot(plotdf, aes(x=cluster,y=variable,fill=value))+
        geom_tile(na.rm=F)+
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(),
              axis.text.x  = element_text(colour = "black", size=12),
              axis.text.y  = element_text(colour = "black", size=12),
              axis.title.x = element_text(colour = "black", size=12, margin = margin(10,0,0,0)), 
              axis.title.y = element_blank(),
              legend.title = element_text(colour = "black", size=12)) +
        labs(x = "cluster N", y = "label") +
        scale_fill_manual(drop=FALSE, values=colorRampPalette(c("yellow","red"))(6), na.value="white", name="count") +
        coord_fixed())
dev.off()

data_card <- matrix(data=NA,nrow=10,ncol=10)
for (clusterA in 1:10){
  for (clusterB in 1:10){
    data_card[clusterA,clusterB] = sum(Ldata[,clusterA+3]*Ldata[,clusterB+3])
  }
}
for(x in 1:10){for(y in 1:10){ if (x<y){data_card[x,y]=NA}}} # remove duplicate values
df = data.frame(cluster = clusters, t(data_card))
colnames(df)[2:11] <- clusters
plotdf = melt(df, id = "cluster")
plotdf$value = cut(plotdf$value, breaks=27, labels=0:26, include.lowest = T)
plotdf$value <- factor(plotdf$value, levels = rev(levels(plotdf$value)))
plotdf$variable <- factor(plotdf$variable, levels = rev(levels(plotdf$variable)))
pdf(file = file.path(projDir, "report/figures/wordle/matrix_dp.pdf"),
    bg="transparent",
    width=5, 
    height=4)
print(ggplot(plotdf, aes(x=cluster,y=variable,fill=value))+
        geom_tile(na.rm=F)+
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(),
              axis.text.x  = element_text(colour = "black", size=14),
              axis.text.y  = element_text(colour = "black", size=14),
              axis.title.x = element_text(colour = "black", size=14, margin = margin(10,0,0,0)), 
              axis.title.y = element_text(colour = "black", size=14, margin = margin(0,10,0,0)),
              legend.title = element_text(colour = "black", size=14)) +
        labs(x = "cluster N", y = "cluster N") +
        scale_fill_manual(drop=FALSE, values=colorRampPalette(c("yellow","red"))(27), na.value="white") +
        coord_fixed())
dev.off()

# neural v GIST & card (original method) --------------------------------
# Wang regions
gist_mat <- vector_gist[between_vals] - min(vector_gist[between_vals]); gist_mat = gist_mat/max(gist_mat)
card_mat <- vector_card[between_vals] - min(vector_card[between_vals]); card_mat = card_mat/max(card_mat)

resdir <- file.path(projDir, 'results/MVPA/between/Wang_2015/all_regions')
regressionFile <- file.path(resdir, 'regression.txt')
sink(regressionFile, append = F)
cat('\n### multiple regression of model matrices v neural data ###\n')
sink()
for (cond in conditions){
  cols2use <- get(sprintf('%s_cols', cond))
  plot_df <- data.frame(region = rep(wang_regions, 2), stream = c(rep('posterior',7), rep('ventral',4), rep('lateral',5)), model = rep(c('gist','card'), each = 20), mean = NA, se = NA, t = NA, df = 20, CI_lo = NA, CI_hi = NA, p = NA, p.sig = NA, p.fdr = NA, p.fdr.sig = NA, cohensD = NA)
  comp_df <- data.frame(region = wang_regions, stream = c(rep('posterior',7), rep('ventral',4), rep('lateral',5)), mean = NA, se = NA, t = NA, df = 20, CI_lo = NA, CI_hi = NA, p = NA, p.sig = NA, p.fdr = NA, p.fdr.sig = NA, cohensD = NA)
  for (region in wang_regions){
    gist_betas <- c()
    card_betas <- c()
    for (subj in levels(data$subject)){
      subj_mat <- as.numeric(data[data$subject == subj & data$region == region, cols2use[between_vals]+3])
      Ldata <- data.frame(neural = subj_mat, gist = gist_mat, card = card_mat)
      Lmodel <- lm(neural ~ gist + card, data = Ldata)
      gist_betas <- append(gist_betas, coef(summary(Lmodel))["gist","Estimate"])
      card_betas <- append(card_betas, coef(summary(Lmodel))["card","Estimate"])
    }
    for (model in c('gist','card')){
      betas <- get(sprintf('%s_betas', model))
      plot_df$mean[plot_df$region == region & plot_df$model == model] = mean(betas)
      plot_df$se[plot_df$region == region & plot_df$model == model] = std.error(betas)
      model_test = t.test(betas,mu = 0)
      plot_df$t[plot_df$region == region & plot_df$model == model] = model_test$statistic
      plot_df$CI_lo[plot_df$region == region & plot_df$model == model] = model_test$conf.int[1]
      plot_df$CI_hi[plot_df$region == region & plot_df$model == model] = model_test$conf.int[2]
      plot_df$p[plot_df$region == region & plot_df$model == model] = model_test$p.value
      plot_df$cohensD[plot_df$region == region & plot_df$model == model] = cohensD(betas, mu = 0)
    }
    comp_test = t.test(card_betas, gist_betas, paired = T)
    comp_df$mean[comp_df$region == region] = mean(card_betas-gist_betas)
    comp_df$se[comp_df$region == region] = std.error(card_betas-gist_betas)
    comp_df$t[comp_df$region == region] = comp_test$statistic
    comp_df$CI_lo[comp_df$region == region] = comp_test$conf.int[1]
    comp_df$CI_hi[comp_df$region == region] = comp_test$conf.int[2]
    comp_df$p[comp_df$region == region] = comp_test$p.value
    comp_df$cohensD[comp_df$region == region] = cohensD(card_betas, gist_betas, method = 'paired')
  }
  plot_df$p.sig <- sig_code(plot_df$p)
  plot_df$p.fdr <- p.adjust(plot_df$p, method = 'fdr')
  plot_df$p.fdr.sig <- sig_code(plot_df$p.fdr)
  plot_df$region = factor(plot_df$region, levels = wang_regions)
  plot_df$stream = factor(plot_df$stream, levels = c('posterior', 'ventral', 'lateral','dorsal'))
  plot_df$model = factor(plot_df$model, levels = c('gist','card'))
  
  comp_df$p.sig <- sig_code(comp_df$p)
  comp_df$p.fdr <- p.adjust(comp_df$p, method = 'fdr')
  comp_df$p.fdr.sig <- sig_code(comp_df$p.fdr)

  sink(regressionFile, append = T)
  cat(sprintf('\n### %s ###\n', cond))
  cat(capture.output(plot_df), sep ='\n')
  cat('\n### card v gist ###\n')
  cat(capture.output(comp_df), sep ='\n')
  sink()
  
  # plot separately for each model
  for (model in c('gist','card')){
    if (model == 'gist'){ ylims = c(-.07, .65)}
    if (model == 'card'){ ylims = c(-.1, .4)}
    temp_plot_df <- plot_df[plot_df$model == model,]
    pdf(file=sprintf('%s/multReg_%s_%s_performance.pdf', resdir, cond, model),
        bg="transparent",
        width=9, 
        height=2.5)
    print(ggplot(temp_plot_df, aes(x=region, y=mean)) + 
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(), 
                  panel.background = element_blank(), 
                  axis.line.x = element_line(colour = "black", size = .5), 
                  axis.line.y = element_line(colour = "black", size = .5), 
                  axis.ticks = element_line(colour = "black", size = .5),
                  axis.text.x  = element_text(colour = "black", size=10), 
                  axis.text.y  = element_text(colour = "black", size=10), 
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(size=13, margin = margin(0,10,0,0)),
                  strip.text = element_text(colour = 'black', size = 13),
                  legend.text = element_text(colour = "black", size=13),
                  legend.title = element_blank()) +
            geom_errorbar(aes(ymin=mean-se, ymax=mean+se), stat = "identity", width=.3, size = .5) +
            geom_bar(stat = "identity", position = 'dodge', width = .7, colour = "black", fill = 'white') +
            facet_grid(.~stream, scales="free", space="free_x", shrink=T) +
            geom_text(aes(label = p.sig, y = mean+se+.02), colour = "black", size = 5) +
            geom_hline(yintercept = 0) +
            scale_y_continuous(limits = ylims, breaks = seq(0,.6,.2)) +
            labs(y = expression(paste("regression coefficient (", beta, ")", sep = ""))))
    dev.off()
  }
}
  
# ventral stream
regressionFile <- file.path(projDir, 'results/MVPA/between/HOCSA/ventral_stream/regression.txt')
sink(regressionFile, append = F)
cat('### multiple regression of model matrices v neural data ###\n')
sink()
for (cond in conditions){
  cols2use <- get(sprintf('%s_cols', cond))
  plot_df <- data.frame(model = c('gist','card'), mean = NA, se = NA, t = NA, df = 20, CI_lo = NA, CI_hi = NA, p = NA, p.sig = NA, p.fdr = NA, p.fdr.sig = NA, cohensD = NA)
  gist_betas <- c()
  card_betas <- c()
  for (subj in levels(data$subject)){
    subj_mat <- as.numeric(data[data$subject == subj & data$region == 'ventral_stream', cols2use[between_vals]+3])
    Ldata <- data.frame(neural = subj_mat, gist = vector_gist, card = as.numeric(na.omit(vector_card[between_vals])))
    Lmodel <- lm(neural ~ gist + card, data = Ldata)
    gist_betas <- append(gist_betas, coef(summary(Lmodel))["gist","Estimate"])
    card_betas <- append(card_betas, coef(summary(Lmodel))["card","Estimate"])
  }
  for (model in c('gist','card')){
    betas <- get(sprintf('%s_betas', model))
    plot_df$mean[plot_df$model == model] = mean(betas)
    plot_df$se[plot_df$model == model] = std.error(betas)
    model_test = t.test(betas,mu = 0)
    plot_df$t[plot_df$model == model] = model_test$statistic
    plot_df$CI_lo[plot_df$model == model] = model_test$conf.int[1]
    plot_df$CI_hi[plot_df$model == model] = model_test$conf.int[2]
    plot_df$p[plot_df$model == model] = model_test$p.value
    plot_df$cohensD[plot_df$model == model] = cohensD(betas, mu = 0)
  }
  plot_df$p.sig <- sig_code(plot_df$p)
  plot_df$p.fdr <- p.adjust(plot_df$p, method = 'fdr')
  plot_df$p.fdr.sig <- sig_code(plot_df$p.fdr)
  
  sink(regressionFile, append = T)
  cat(sprintf('\n### %s ###\n', cond))
  cat(capture.output(plot_df), sep ='\n')
  sink()
  
  # compare GIST and card
  model_test <- t.test(gist_betas, card_betas, paired = T)
  model_d = cohensD(gist_betas, card_betas, method = 'paired')
  sink(regressionFile, append = T)
  cat('\n### gist v card ###\n')
  cat(capture.output(model_test), sep ='\n')
  cat("Cohen's d")
  cat(capture.output(model_d), sep ='\n')
  sink()
  
  # run partial correlations between neural and GIST controlling for a categorical model
  catModel <- c(1,0,1,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1)
  gistModel <- vector_gistWB
  GISTctrlSem <- c()
  semCtrlGIST <- c()
  for (subj in levels(data$subject)){
    subj_mat <- as.numeric(na.omit(as.numeric(data[data$subject == subj & data$region == 'ventral_stream', cols2use+3])))
    GISTctrlSemTest <- pcor.test(subj_mat, gistModel, catModel)
    semCtrlGISTtest <- pcor.test(subj_mat, catModel, gistModel)
    GISTctrlSem <- c(GISTctrlSem, GISTctrlSemTest$estimate)
    semCtrlGIST <- c(semCtrlGIST, semCtrlGISTtest$estimate)
  }
  
  sink(regressionFile, append = T)
  cat('\n### Tim-requested partial correlation ###\n')
  cat('\n## GIST controlling for W v B ##\n')
  cat(capture.output(t.test(GISTctrlSem, mu = 0)), sep='\n')
  cat('\n## W v B controlling for GIST ##\n')
  cat(capture.output(t.test(semCtrlGIST, mu = 0)), sep='\n')
  sink()
  
  # do step-wise regression version of previous analysis
  GISTctrlSem <- c()
  semCtrlGIST <- c()
  GISTnorm <- c()
  semNorm <- c()
  for (subj in levels(data$subject)){
    subj_mat <- as.numeric(na.omit(as.numeric(data[data$subject == subj & data$region == 'ventral_stream', cols2use+3])))
    
    # semantic after removing GIST
    model = lm(subj_mat ~ gistModel)
    semNormTest = summary(model)
    semNorm <- c(semNorm, semNormTest$coefficients[2])
    residual = residuals(model)
    semCtrlGISTtest = summary(lm(residual ~ catModel))
    semCtrlGIST <- c(semCtrlGIST, semCtrlGISTtest$coefficients[2])
    
    # GIST after removing semantic
    model = lm(subj_mat ~ catModel)
    GISTnormTest = summary(model)
    GISTnorm <- c(GISTnorm, GISTnormTest$coefficients[2])
    residual = residuals(model)
    GISTctrlSemTest = summary(lm(residual ~ gistModel))
    GISTctrlSem <- c(GISTctrlSem, GISTctrlSemTest$coefficients[2])
  }
  
  sink(regressionFile, append = T)
  cat('\n### Tim-requested step-wise regression ###\n')
  cat('\n## GIST ##\n')
  cat(capture.output(t.test(GISTnorm, mu = 0)), sep='\n')
  cat('\n## W v B ##\n')
  cat(capture.output(t.test(semNorm, mu = 0)), sep='\n')
  cat('\n## GIST controlling for W v B ##\n')
  cat(capture.output(t.test(GISTctrlSem, mu = 0)), sep='\n')
  cat('\n## W v B controlling for GIST ##\n')
  cat(capture.output(t.test(semCtrlGIST, mu = 0)), sep='\n')
  sink()
}


# neural v GIST & card (permutation method) --------------------------------
# NOTE: for this analysis, ventral stats are stored with wang regions in text file stated below.
regressionFile <- file.path(resdir, 'regression_perm.txt')
sink(regressionFile, append = F)
cat('### multiple regression of model matrices v neural data ###\n')
sink()

gistMatrix <- vector_gist[between_vals] - min(vector_gist[between_vals]); gist_mat = gist_mat/max(gist_mat)
cardMatrix <- vector_card[between_vals] - min(vector_card[between_vals]); card_mat = card_mat/max(card_mat)

for (cond in conditions){
  cols2use <- get(sprintf('%s_cols', cond))
  plot_df <- data.frame(region = rep(wang_regions, 2), stream = c(rep('posterior',7), rep('ventral',4), rep('lateral',5)), model = rep(c('gist','card'), each = 20), mean = NA, se = NA, t = NA, df = 20, CI_lo = NA, CI_hi = NA, p = NA, p.sig = NA, p.fdr = NA, p.fdr.sig = NA, cohensD = NA)
  comp_df <- data.frame(region = wang_regions, stream = c(rep('posterior',7), rep('ventral',4), rep('lateral',5)), mean = NA, se = NA, t = NA, df = 20, CI_lo = NA, CI_hi = NA, p = NA, p.sig = NA, p.fdr = NA, p.fdr.sig = NA, cohensD = NA)
  for (region in wang_regions){
    gistBetas <- c()
    cardBetas <- c()
    for (subj in levels(data$subject)){
      shuffledGistBetas <- c()
      shuffledCardBetas <- c()
      neuralMatrix <- as.numeric(data[data$subject == subj & data$region == region, cols2use[between_vals]+3])
      for (p in 1:1000){
        shuffledMatrix <- sample(neuralMatrix, size = length(neuralMatrix), replace = F)
        df <- data.frame(neural = shuffledMatrix, gist = gistMatrix, card = cardMatrix)
        model <- lm(neural ~ gist + card, data = df)
        shuffledGistBetas <- append(shuffledGistBetas, coef(summary(model))["gist","Estimate"])
        shuffledCardBetas <- append(shuffledCardBetas, coef(summary(model))["card","Estimate"])
      }
      df <- data.frame(neural = neuralMatrix, gist = gistMatrix, card = cardMatrix)
      model <- lm(neural ~ gist + card, data = df)
      originalGistBeta <- coef(summary(model))["gist","Estimate"]
      originalCardBeta <- coef(summary(model))["card","Estimate"]
      gistBetas <- append(gistBetas, originalGistBeta - mean(shuffledGistBetas))
      cardBetas <- append(cardBetas, originalCardBeta - mean(shuffledCardBetas))
    }
    for (model in c('gist','card')){
      betas <- get(sprintf('%sBetas', model))
      plot_df$mean[plot_df$region == region & plot_df$model == model] = mean(betas)
      plot_df$se[plot_df$region == region & plot_df$model == model] = std.error(betas)
      model_test = t.test(betas,mu = 0)
      plot_df$t[plot_df$region == region & plot_df$model == model] = model_test$statistic
      plot_df$CI_lo[plot_df$region == region & plot_df$model == model] = model_test$conf.int[1]
      plot_df$CI_hi[plot_df$region == region & plot_df$model == model] = model_test$conf.int[2]
      plot_df$p[plot_df$region == region & plot_df$model == model] = model_test$p.value
      plot_df$cohensD[plot_df$region == region & plot_df$model == model] = cohensD(betas, mu = 0)
    }
    comp_test = t.test(cardBetas, gistBetas, paired = T)
    comp_df$mean[comp_df$region == region] = mean(cardBetas-gistBetas)
    comp_df$se[comp_df$region == region] = std.error(cardBetas-gistBetas)
    comp_df$t[comp_df$region == region] = comp_test$statistic
    comp_df$CI_lo[comp_df$region == region] = comp_test$conf.int[1]
    comp_df$CI_hi[comp_df$region == region] = comp_test$conf.int[2]
    comp_df$p[comp_df$region == region] = comp_test$p.value
    comp_df$cohensD[comp_df$region == region] = cohensD(cardBetas, gistBetas, method = 'paired')
  }
  plot_df$p.sig <- sig_code(plot_df$p)
  plot_df$p.fdr <- p.adjust(plot_df$p, method = 'fdr')
  plot_df$p.fdr.sig <- sig_code(plot_df$p.fdr)
  plot_df$region = factor(plot_df$region, levels = wang_regions)
  plot_df$stream = factor(plot_df$stream, levels = c('posterior', 'ventral', 'lateral','dorsal'))
  plot_df$model = factor(plot_df$model, levels = c('gist','card'))
  
  comp_df$p.sig <- sig_code(comp_df$p)
  comp_df$p.fdr <- p.adjust(comp_df$p, method = 'fdr')
  comp_df$p.fdr.sig <- sig_code(comp_df$p.fdr)
  
  sink(regressionFile, append = T)
  cat(sprintf('\n### %s ###\n', cond))
  cat(capture.output(plot_df), sep ='\n')
  cat('\n### card v gist ###\n')
  cat(capture.output(comp_df), sep ='\n')
  sink()
}


# MDS for V1, ventral, card and gist --------------------------------------

library(plotrix)
# get matrices
gist <- matrix_gist
card <- matrix_card
V1 <- matrix(data = as.numeric(colMeans(data[data$region == 'V1v', grey_cols+3])), nrow = 10, byrow = T)
ventral <- matrix(data = as.numeric(colMeans(data[data$region == 'ventral_stream', grey_cols+3])), nrow = 10, byrow = T)

matricesMDS <- c('gist', 'card', 'V1','ventral')
for (m in 1:4){
  ma <- get(matricesMDS[m])
  for (r in 1:10){
    for (c in 1:10){
      ma[r,c] <- ma[c,r]
    }
  }
  ma = plotrix::rescale(ma, c(0,1))
  ma = 1-ma
  row.names(ma) = clusters
  colnames(ma) = clusters
  fit = cmdscale(ma, eig = T, k = 2)
  x = fit$points[,1]
  y = fit$points[,2]
  mds_plot_data = data.frame(x = x, y = y)
  pdf(file = paste(projDir, "/report/figures/dump/MDS_", matricesMDS[m], ".pdf", sep=""),
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











# region matrix correlation across exp1 and exp2 --------------------------
cors_vect2 <- read.csv('/Volumes/ddc/research/projects/p007/region_matrix_exp2.csv')
inv_cors_mat.avdiag = inv_cors_mat
for (r in 1:16){ for (c in 2:17){ if (r+1 <= c) inv_cors_mat.avdiag[r,c] <- NA}}
cors_vect <- c(na.omit(c(as.matrix(inv_cors_mat.avdiag[,-1]))))
write.csv(cors_vect, '/Volumes/ddc/research/projects/p007/region_matrix_exp1.csv')

statsFile <- '/Volumes/ddc/research/projects/p007/region_matrix_stats.txt'
sink(statsFile)
cat('### correlation between off-diagonal elements for region matrices in experiments 1 and 2 ###\n')
cat(capture.output(cor.test(cors_vect, cors_vect2$x)), sep = '\n')
sink()

