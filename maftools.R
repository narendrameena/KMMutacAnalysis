#author narumeena 
#description analysis of vcf file of cpc 2 using maftools 

suppressPackageStartupMessages(require(maftools))


#read annovar file into maf formate 
maf = annovarToMaf(annovar = "data/mutak2.GATK.mutac.raw.vcf.hg19_multianno.txt", Center = 'BISR', refBuild = 'hg19', 
                               tsbCol = 'Tumor_Sample_Barcode', table = 'ensGene')



laml = read.maf(maf = maf, removeSilent = TRUE, useAll = FALSE)

#Typing laml shows basic summary of MAF file.
laml

#Shows sample summry.
getSampleSummary(laml)


#Shows frequently mutated genes.
getGeneSummary(laml)


#Shows all fields in MAF
getFields(laml)

#Writes maf summary to an output file with basename laml.
write.mafSummary(maf = laml, basename = 'data/laml')

plotmafSummary(maf = laml, rmOutlier = FALSE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE,top=35,height=1000)


#We will draw oncoplots for top ten mutated genes. (Removing non-mutated samples from the plot for better visualization)
##example 
#laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#laml <- read.maf(maf = laml.maf, removeSilent = TRUE, useAll = FALSE)
#oncoplot(maf = laml, top = 35)


oncoplot(maf = laml, top = 35, removeNonMutated = TRUE)


#read TCGA maf file for LAML
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
all.lesions <- system.file("extdata", "all_lesions.conf_99.txt", package = "maftools")
amp.genes <- system.file("extdata", "amp_genes.conf_99.txt", package = "maftools")
del.genes <- system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
laml.plus.gistic = read.maf(maf = laml.maf, removeSilent = TRUE, useAll = FALSE, gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, isTCGA = TRUE)


#We will draw oncoplots for top ten mutated genes. (Removing non-mutated samples from the plot for better visualization)
oncoplot(maf = laml.plus.gistic, top = 10, removeNonMutated = TRUE, sortByMutation = TRUE)


#Read FAB classification of TCGA LAML barcodes.
laml.fab.anno = system.file('extdata', 'tcga_laml_fab_annotation.txt', package = 'maftools')
laml.fab.anno = read.delim(laml.fab.anno, sep = '\t')
head(laml.fab.anno)


#Changing colors (You can use any colors, here in this example we will use a color palette from RColorBrewer)
col = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')
#We will plot same top ten mutated genes with FAB classification as annotation and using above defined colors.
oncoplot(maf = laml, top = 10, annotation = laml.fab.anno, removeNonMutated = TRUE, colors = col)




oncostrip(maf = laml, genes = c('DNMT3A','NPM1', 'RUNX1'), removeNonMutated = TRUE, showTumorSampleBarcodes = FALSE)


laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)


#Lets plot lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
dnmt3a.lpop = lollipopPlot(maf = laml, gene = 'TTN', AACol = 'Protein_Change', showMutationRate = TRUE, domainLabelSize = 3, defaultYaxis = FALSE)
#


#Lets plot mutations on KIT gene, without repel option.
kit.lpop = lollipopPlot(maf = laml, gene = 'TTN', AACol = 'Protein_Change', labelPos = c(416, 418), refSeqID = 'NM_000222', domainLabelSize = 3)


#Same plot with repel=TRUE
kit.lpop = lollipopPlot(maf = laml, gene = 'TTN', AACol = 'Protein_Change', labelPos = c(416, 418), refSeqID = 'NM_000222', repel = TRUE, domainLabelSize = 3)


laml.dnmt3a = lollipopPlot(maf = laml, gene = 'TTN', AACol = 'Protein_Change', refSeqID = 'NM_175629', labelPos = 882, collapsePosLabel = TRUE, cBioPortal = TRUE, domainLabelSize = 3, defaultYaxis = FALSE)


tcga.ab.009.seg <- system.file("extdata", "cpc2", package = "maftools")

           verbose = TRUE
plotCBSsegments(cbsFile = tcga.ab.009.seg, maf = laml, labelAll = TRUE)


coad <- system.file("extdata", "coad.maf.gz", package = "maftools")
coad = read.maf(maf = coad)

coad.rf = rainfallPlot(maf = laml, detectChangePoints = TRUE, fontSize = 12, pointSize = 0.6)


laml.mutload = tcgaCompare(maf = laml, cohortName = 'cpc2')


geneCloud(input = laml, minMut = 3)

all.lesions <- system.file("extdata", "all_lesions.conf_99.txt", package = "maftools")
amp.genes <- system.file("extdata", "amp_genes.conf_99.txt", package = "maftools")
del.genes <- system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
laml.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, isTCGA = TRUE)
#


gisticPlot(gistic = laml.gistic)

plotGisticResults(gistic = laml.gistic)


#We will run mutExclusive on top 10 mutated genes. 
laml.mut.excl = mutExclusive(maf = laml, top = 1)
head(laml.mut.excl)


oncostrip(maf = laml, genes = c('NPM1', 'RUNX1'), sort = TRUE, removeNonMutated = TRUE)

laml.sig = oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
head(laml.sig)


plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE)


laml.pfam = pfamDomains(maf = laml, AACol = 'Protein_Change', top = 3)


#Protein summary (Printing first 7 columns for display convenience)
laml.pfam$proteinSummary[,1:7, with = FALSE]


#Domain summary (Printing first 3 columns for display convenience)
laml.pfam$domainSummary[,1:3, with = FALSE]

#MutsigCV results for TCGA-AML
laml.mutsig <- system.file("extdata", "LAML_sig_genes.txt.gz", package = "maftools")
laml.pancan = pancanComparision(mutsigResults = laml.mutsig, qval = 0.1, cohortName = 'LAML', inputSampleSize = 200, label = 1, normSampleSize = TRUE)
#


laml.surv <- system.file("extdata", "laml_survival.tsv", package = "maftools")
laml.surv = read.delim(file = laml.surv, header = TRUE, stringsAsFactors = FALSE)
head(laml.surv)


#Survival analysis based on grouping of DNMT3A mutation status
laml.dnmt3a.surv = mafSurvival(maf = laml, clinicalData = laml.surv, genes = 'TTN', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', showConfInt = TRUE)
#
