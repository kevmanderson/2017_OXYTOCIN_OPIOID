# Supplemental code for 
# 
# Boosting the social effects of oxytocin with opioid antagonism 
#
# Olga Dal Monte*1, Matthew Piva*1,2,3, Kevin Anderson1, Marios Tringides1, Avram Holmes1,4, Steve W. C. Chang1,2,3
# 1 Department of Psychology, Yale University, New Haven, CT 06520
# 2 Department of Neuroscience, Yale University School of Medicine, New Haven, CT 06510
# 3 Interdepartmental Neuroscience Program, Yale University School of Medicine, New Haven, CT 06510
# 4 Department of Psychiatry, Yale University School of Medicine, New Haven, CT 06510
#
# Contact: Kevin M. Anderson, kevin.anderson@yale.edu 


# SET BASE DIRECTORY; source function library
# ----------------
base.dir     <- '/Users/kevinanderson/PHD/PROJECTS/OXYTOCIN/'
function.lib <- paste(base.dir, 'scripts/oxy_function_library.R', sep = '')
source(function.lib)

# LOAD PACKAGES
# ----------------
packages.to.load <- paste(base.dir, 'scripts/packages_to_load.txt', sep = '')
loadPackages(packages.to.load)
# ----------------


# SUBJECT LIST
# ----------------
donor.nums <- c('9861', '10021', '12876', '14380', '15496', '15697')
# ----------------


# PRE-PROCESS MICROARRAY DATA
# ----------------
data_path <- paste(base.dir, 'lib/AHBA', sep = '')

# TRUE/FALSE; either pre-process the data from scratch (takes a few minutes), or
# read in previously pre-processed saved Rdata object
preprocess <- FALSE; #TRUE

# DON'T PREPROCESS
if (preprocess == FALSE){
  load_path <- paste(data_path, '/all_data.Rdata', sep = '')
  load(load_path)
  
# PREPROCESS
} else {

  # Read/process data from each subject
  # -------------------------------------------
  all_data <- NULL
  for ( donor in donor.nums ) {
    file <- paste('donor_', donor, sep='')
    print(paste('Reading and collapsing data for: ', donor, sep = ''))
    
    all_data[[donor]] <- NULL
    
    # Sample Information
    # ---------------fil----
    saname    <- paste(data_path, file, 'SampleAnnot.csv', sep='/')
    samp_info <- read.csv(saname)
    all_data[[donor]]$raw_samps   <- samp_info 
    
    # Read Microexpression Data
    # -------------------------------------------------------
    fname                <- paste(data_path, file, 'MicroarrayExpression.csv', sep='/')
    microdata            <- fread(fname, header = F, sep = ',')
    micro_arr            <- as.matrix(microdata)
    
    micro_temp           <- micro_arr[,2:dim(micro_arr)[2]] # First column contains probe IDs
    rownames(micro_temp) <- micro_arr[,1]
    micro_df             <- as.data.frame(micro_temp)
    
    all_data[[donor]]$raw_micros <- micro_df

    # Read PA-Call (signal present vs absent)
    # -------------------------------------------------------
    paname      <- paste(data_path, file, 'PACall.csv', sep='/')
    pacall      <- fread(paname, header = F, sep = ',')
    pacall_arr  <- as.matrix(pacall)
    pa_dat      <- pacall_arr[,2:dim(pacall_arr)[2]]
    rownames(pa_dat) <- pacall_arr[,1]
    all_data[[donor]]$raw_pas <- pa_dat
    
    # Information about each Gene Probe (~50,000 probes for ~20,000 genes)
    # -------------------------------------------------------
    pname     <- paste(data_path, file, 'Probes.csv', sep='/')
    probes    <- read.csv(pname)
    all_data[[donor]]$raw_probes <- probes
    
    # Gene Ontology Info
    # -------------------------------------------------------
    oname    <- paste(data_path, file, 'Ontology.csv', sep='/')
    ont_data <- read.csv(oname)
    all_data[[donor]]$raw_ont <- ont_data
    
    # Get rid of probes without an entrez-id
    # -------------------------------------------------------
    num_samples   <- dim(all_data[[donor]]$raw_micros)[2]
    trash_me      <- is.na(all_data[[donor]]$raw_probes$entrez_id) # & pasums > cutoff
    good_probes   <- all_data[[donor]]$raw_probes[trash_me == FALSE,]  
    
    # Select the good probes 
    # -------------------------------------------------------
    all_data[[donor]]$probes_filter <- good_probes
    all_data[[donor]]$micro_filter  <- all_data[[donor]]$raw_micros[rownames(all_data[[donor]]$raw_micros) %in% good_probes$probe_id,]
    all_data[[donor]]$pas_filter    <- all_data[[donor]]$raw_pas[rownames(all_data[[donor]]$raw_pas) %in% good_probes$probe_id,]
    
    # Collapse duplicate probes
    # -------------------------------------------------------
    agilent.select <- read.csv( paste(base.dir, 'lib/AgilentProbeSelection.csv', sep = '') )
    agilent.filter <- agilent.select[!is.na(agilent.select$entrez_id),]
    use.probe      <- agilent.filter$probe_name[agilent.filter$is_best_probe]
    all_data[[donor]]$micro_filter[which(all_data[[donor]]$probes_filter$probe_name %in% use.probe),]

    
    # Collect the output
    # -------------------------------------------------------
    all_data[[donor]]$micro_collapsed  <- all_data[[donor]]$micro_filter[which(all_data[[donor]]$probes_filter$probe_name %in% use.probe),]
    all_data[[donor]]$probes_collapse  <- all_data[[donor]]$probes_filter[which(all_data[[donor]]$probes_filter$probe_name %in% use.probe),]
    all_data[[donor]]$pas_collapse     <- all_data[[donor]]$pas_filter[which(all_data[[donor]]$probes_filter$probe_name %in% use.probe),]
    
  }
  # Save the pre-processed data for later loading
  # -------------------------------------------------------
  save_path <- paste(data_path, '/all_data.Rdata', sep = '')
  save(file = save_path, x = all_data)
}


# Calculate information about the frequency of regions across each of the six donors
# ----------------
sample.array <- NULL
name.array   <- NULL
for (donor in donor.nums){
  rownames(all_data[[donor]]$micro_collapsed) <- all_data[[donor]]$probes_collapse$gene_symbol
  
  cur.samples    <- all_data[[donor]]$raw_samps$structure_acronym #sample acronyms that are present for this subject
  unique.samples <- as.character(unique(cur.samples)) 
  sample.array   <- c(sample.array, unique.samples) 
  cur.names      <- all_data[[donor]]$raw_samps$structure_name 
  unique.names   <- as.character(unique(gsub('left|right','', cur.names))) # strip left/right
  name.array     <- c(name.array, unique.names)
}
region.frequency <- rev(sort(table(sample.array))) 
# Only examine regions with samples from at least 4 donors
use.regions      <- names(region.frequency[region.frequency >= 4]) 



# Ontology information regarding the structure hieararchy of the sample
# -----------------
ontology <- all_data$`9861`$raw_ont 


out.names <- NULL
for (reg in use.regions){
  cur.name <- ontology$name[which(ontology$acronym == reg)[1]]
  out.name <- gsub(', left|, right', '', as.character(cur.name))
  out.names <- c(out.names, out.name)
}


# For each subject, average duplicate samples in each of the 190 regions
# ----------------
ct   <- 1
all.regional.gene.expr <- NULL
for ( donor in donor.nums ){ # for each donor 
  print(donor)
  regional.gene.expr <- NULL
  col.acros          <- all_data[[donor]]$raw_samps$structure_acronym

  # Each region defined above
  for (acro in use.regions){ 
    acro.idxs  <- which(col.acros == acro) # get indexes of matching acronyms for subject and acronyms where frequency is at least 4
    if (length(acro.idxs) == 0){ # not present
      reg.expr <- NA 
    } else if (length(acro.idxs) > 1){
      reg.expr <- rowMeans(all_data[[donor]]$micro_collapsed[, acro.idxs]) 
    } else {
      reg.expr <- all_data[[donor]]$micro_collapsed[,acro.idxs] # if one sample, use genes in that region for this subject
    }
    regional.gene.expr <- cbind(regional.gene.expr, reg.expr) # put data in regional.gene.expr
  }
  colnames(regional.gene.expr) <- use.regions 
  all.regional.gene.expr[[ct]] <- regional.gene.expr
  ct <- ct + 1
}


# Average across subjects
# ----------------
avg.expr.matrix <- NULL
for ( acro in use.regions ){ 
  print(acro)
  col.idx <- which(colnames(all.regional.gene.expr[[1]]) == acro) 
  col.dat <- cbind(all.regional.gene.expr[[1]][,col.idx], all.regional.gene.expr[[4]][,col.idx], 
                   all.regional.gene.expr[[2]][,col.idx], all.regional.gene.expr[[5]][,col.idx],  
                   all.regional.gene.expr[[3]][,col.idx], all.regional.gene.expr[[6]][,col.idx]) 
  avg.expr.matrix <- cbind(avg.expr.matrix, rowMeans(col.dat, na.rm = TRUE)) # average data across subjects
}
colnames(avg.expr.matrix) <- use.regions # name columns with acronyms
# ----------------


# Get data ready for plotting
# ----------------
plot.order.path <- paste(base.dir, 'data/plot_order.csv', sep = '')
plot.order      <- read.csv(plot.order.path, header = FALSE) # create data frame from file


# Get data ready for plotting
# ----------------
gene             <- 'OXT' # 'OXTR', 'AVP', 'AVPR1A', 'OPRM1', 'OPRK1', 'OPRD1'
plot.expr        <- avg.expr.matrix[rownames(avg.expr.matrix) == gene] # find OXTR
names(plot.expr) <- colnames(avg.expr.matrix) # name plot.expr with acronyms from use.regions 
ordered.regional.gene.expr <- NULL
color.arr        <- NULL
name.arr         <- NULL
ontology.splits  <- strsplit2(ontology$structure_id_path, '/') # split acronyms in ontology


# Put regions in order of their overarching ontology
# ----------------
unique.colors <- NULL
ct <- 1
name.arr <- NULL
color.arr <- NULL
ordered.regional.gene.expr <- NULL
for (ontology.number in plot.order$V2){ # for each ontology number
  plot.gene.expr.mean <- plot.expr
  match.expr          <- paste('\\b', ontology.number, '\\b', sep = '') # add \\b to front and back of ontology number
  
  ont.matches <- NULL
  for (row in 1:dim(ontology.splits)[1] ){ # for each row until the number of acronyms in ontology.splits
    if (length(grep(match.expr, ontology.splits[row,])) == 0){ # if not target ontology number  
      ont.matches <- c(ont.matches, 0) # set ont.matches to 0
    } else {
      ont.matches <- c(ont.matches, grep(match.expr, ontology.splits[row,])) # if target ontology number, get index of matching ontology number
    }
  }

  acros.in.this.category     <- as.character(ontology$acronym[which(ont.matches > 0)]) 
  row.idxs                   <- which(use.regions %in% acros.in.this.category) 
  if (length(row.idxs) > 0){
    unique.colors <- c(unique.colors, as.character(plot.order$V3[ct]))
  }
  ordered.regional.gene.expr <- c(ordered.regional.gene.expr, plot.expr[row.idxs]) # get these acronyms for target gene
  color.arr <- c(color.arr, as.character(rep(plot.order$V3[ct], length(row.idxs)))) # get colors, repeat for number of matching acronyms
  name.arr  <- c(name.arr, as.character(rep(plot.order$V1[ct], length(row.idxs)))) # get region categories, repeat for number of matching acronyms
  
  ct <- ct + 1
}
# ----------------
out.name <- paste(base.dir, 'data/region_info.csv', sep = '')
tmp <- data.frame(cbind(name.arr, out.names, use.regions, color.arr))
colnames(tmp) <- c('Structure Category','Structure Name','Acronym', 'HEX')
write.csv(file=out.name, x=tmp, row.names=FALSE)


# make data frame with region, expression, color, and category
# ----------------
plot.me           <- as.data.frame(cbind(names(ordered.regional.gene.expr), ordered.regional.gene.expr, color.arr, name.arr))
plot.me$ordered.regional.gene.expr <-  as.numeric(as.character(plot.me$ordered.regional.gene.expr)) # convert expression to number
plot.me$V1        <- as.character(plot.me$V1)
colnames(plot.me) <- c('region', 'expr', 'color', 'category') # name columns of plot.me
plot.me$expr      <- plot.me$expr - mean(plot.me$expr) # mean normalize
plot.me$category  <- factor(x = as.character(plot.me$category), levels = as.character(plot.order$V1)) # converts category to factor, names levels
rownames(plot.me) <- plot.me$region # names rows of plot.me with region names
# ----------------


ggplot(data=plot.me, aes(x=region, y=expr, fill=category )) + #ylim(-4,4) + 
  geom_bar(stat="identity") +# ylim(-4, 4) + 
  scale_x_discrete(limits = plot.me$region) + 
  scale_fill_manual(values=unique.colors) + 
  geom_hline(yintercept=mean(plot.me$expr) + sd(plot.me$expr)) +
  ggtitle(gene) 


# region-wise exprssion for genes of interest
avp.expr    <- avg.expr.matrix[rownames(avg.expr.matrix) == 'AVP']
avp1ra.expr <- avg.expr.matrix[rownames(avg.expr.matrix) == 'AVPR1A']
oxt.expr    <- avg.expr.matrix[rownames(avg.expr.matrix) == 'OXT']
oprm1.expr  <- avg.expr.matrix[rownames(avg.expr.matrix) == 'OPRM1']
oprk1.expr  <- avg.expr.matrix[rownames(avg.expr.matrix) == 'OPRK1']
oprd1.expr  <- avg.expr.matrix[rownames(avg.expr.matrix) == 'OPRD1']
oxtr.expr   <- avg.expr.matrix[rownames(avg.expr.matrix) == 'OXTR']

# regions 1 standard deviation above mean - OXT
scaled <- as.data.frame(scale(oxt.expr))
oxt.top.regions  <- use.regions[scaled >= 1]
cor.test(oxt.expr[scaled >= 1], oprm1.expr[scaled >= 1])
cor.test(oxt.expr[scaled >= 1], oprk1.expr[scaled >= 1])
cor.test(oxt.expr[scaled >= 1], oprd1.expr[scaled >= 1])

# OPRM1
color.arr  <- rep('grey', 1, length(oxt.expr))
color.arr[scaled >= 1] <- 'red'
plot.me <- data.frame(oxt.expr, oprm1.expr, color.arr, colnames(avg.expr.matrix))
ggplot(plot.me) + geom_point(aes(x=oxt.expr, y=oprm1.expr, color=color.arr)) + scale_color_manual(values=c('darkgrey','red')) + theme_bw() + 
  geom_hline(yintercept = mean(oprm1.expr))

# OPRK1
plot.me <- data.frame(oxt.expr, oprk1.expr, color.arr, colnames(avg.expr.matrix))
ggplot(plot.me) + geom_point(aes(x=oxt.expr, y=oprk1.expr, color=color.arr)) + scale_color_manual(values=c('darkgrey','red')) + theme_bw() + 
  geom_hline(yintercept = mean(oprk1.expr))

# OPRD1
plot.me <- data.frame(oxt.expr, oprd1.expr, color.arr, colnames(avg.expr.matrix))
ggplot(plot.me) + geom_point(aes(x=oxt.expr, y=oprd1.expr, color=color.arr)) + scale_color_manual(values=c('darkgrey','red')) + theme_bw() + 
  geom_hline(yintercept = mean(oprd1.expr))

# OXTR
plot.me <- data.frame(oxt.expr, oxtr.expr, color.arr, colnames(avg.expr.matrix))
ggplot(plot.me) + geom_point(aes(x=oxt.expr, y=oxtr.expr, color=color.arr)) + scale_color_manual(values=c('darkgrey','red')) + theme_bw() + 
  geom_hline(yintercept = mean(oxtr.expr))


# Compare relative expression of Mu, Kappa, and Delta in each OXT enriched region
for (oxt.reg in oxt.top.regions){
  print(oxt.reg)
  col.idx = which(colnames(all.regional.gene.expr[[1]]) == oxt.reg) # get index for where acronym is in this subject
  
  col.dat = cbind(all.regional.gene.expr[[1]][,col.idx], all.regional.gene.expr[[4]][,col.idx], 
                  all.regional.gene.expr[[2]][,col.idx], all.regional.gene.expr[[5]][,col.idx],  
                  all.regional.gene.expr[[3]][,col.idx], all.regional.gene.expr[[6]][,col.idx])
  
  oprm1.dat <- col.dat[rownames(col.dat) == 'OPRM1',]
  oprk1.dat <- col.dat[rownames(col.dat) == 'OPRK1',]
  oprd1.dat <- col.dat[rownames(col.dat) == 'OPRD1',]
  
  # Create a data frame suitable for a within-subjects one way ANOVA
  reg.dat     <- c(oprm1.dat, oprk1.dat, oprd1.dat)
  gene.group  <- c(rep(1,1,6), c(rep(2,1,6), rep(3,1,6)))
  gene.names  <- c(rep('Mu',1,6), c(rep('Kappa',1,6), rep('Delta',1,6)))
  donor.names <- c(donor.nums, donor.nums, donor.nums)
  aov.df  <- data.frame(expr=reg.dat, group=gene.group, names=gene.names, donor=donor.nums)
  aov.df  <- aov.df[!is.na(aov.df$expr),]
  aov.out <- aov(expr ~ group + Error(donor/group), aov.df)
  x <- summary(aov.out)
  print(x$`Error: donor:group`)
  print(pairwise.t.test(aov.df$expr, aov.df$names, p.adj = "fdr"))
  cat('\n\n\n')
  
  pdf(paste('~/PHD/PROJECTS/OXYTOCIN/figures/', oxt.reg, '_OpioidReceptors.pdf', sep = ''))
  
  print(ggplot(aov.df, aes(names, expr)) + geom_point() +
          stat_summary(fun.y = "median", geom = "point", pch = "_", size = 25) + ylim(0,8) +
          geom_line(data = aov.df, aes(x = names, y= expr, group = donor), linetype="dotted") + scale_color_grey() + theme_classic() + 
          ggtitle(oxt.reg) + xlab('Opioid receptor subtype') + ylab('Mean-normalized gene expression (log2)') +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
  
  dev.off()
}












