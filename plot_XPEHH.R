
library(rehh)
require(graphics)

palette(c("#117733","#44aa99"))
dev.off()
filein="refGeneHG19.gz"
dat<-read.table(filein,as.is=T,head=T,comment.char="")

find_genes <- function(x) {
	cat("Chr","start-end Position","","Maximum:","Genes",'\n')
	genes <- c()
	for (g in 1:nrow(cr_wg.xpehh)) {
		chr <- levels(factor(cr_wg.xpehh$CHR[g]))
		maximum <- cr_wg.xpehh$MAX_MRK[g]
		min.pos <- as.integer(cr_wg.xpehh$START[g])
		max.pos <- as.integer(cr_wg.xpehh$END[g])
		xx2 <- dat[dat[,"chrom"]==paste("chr",chr,sep="") & dat[,"cdsStart"]<max.pos & dat[,"cdsEnd"] >min.pos,]
		cat(chr,paste0(min.pos,'-',max.pos),maximum,':',xx2$name2,'\n')
		if (g==1) { genes <- append(genes,c(paste0(xx2$chr,":"),xx2$name2)) }
		else { genes <- append(genes,c(paste0("| ",xx2$chr,":"),xx2$name2)) }
	}
	genes <- genes[!duplicated(genes)]
	genes <- genes[genes!="| :"]
	genes <- genes[genes!=":"]
	return(genes)

}


### iHS for pop1
#Populations_List1 = #c("BurkinaFaso_Fulani_Banfora","BurkinaFaso_Fulani_Tindangou","BurkinaFaso_Fulani_Ziniare","Cameroon_Fulani_Tcheboua","Chad_Fulani_Bongor","Chad_Fulani_Linia","Guinea_Fulani","Mali_Fulani_Diafarabe","Mali_Fulani_InnerDelta","Mauritania_Fulani_Assaba","Niger_Fulani_Abalak","Niger_Fulani_Ader","Niger_Fulani_Balatungur","Niger_Fulani_Diffa","Niger_Fulani_Zinder","Senegal_Fulani_Linguere")


#for (p1 in 1:length(Populations_List1)) {
	#pop1=Populations_List1[[p1]]
	
	#wg.ihs.pop1 <- readRDS(file=paste0("rehh_outputs/",pop1,"_ihsdata.rds"))
	
	# First option to calculate the threshold
	#thres <- unname(15*(quantile(wg.ihs.pop1$ihs$LOGPVALUE, prob=0.99, na.rm=T)))
	#cat('\n',pop1,'iHS',thres,'\n')
	# 
	# Second option to predefine the threshold
	#thres <- 5
	#cat('\n',pop1,'iHS\n')
	#
	#cr <- calc_candidate_regions(wg.ihs.pop1, ignore_sign= T, window_size= 1E4, threshold=thres, pval=T)
	#genes <- find_genes(cr)
	#paste(genes,collapse=' ')
	
	# pdf(file= paste0("rehh_Plots/Edited_iHS_plots/",pop1,"_iHS.plot.pdf"), height=30, width=40, units='cm', res=300)
	#png(paste0("rehh_Plots/Edited_iHS_plots/",pop1,"_iHS.plot.png"), height=30, width=50, units='cm', res=500)
	#par(mfrow=c(2,1))
	#manhattanplot(wg.ihs.pop1, main=paste(genes, collapse=" "), pval= T, pch= 20, threshold= c(-thres,thres), cex= 0.6, adj= 0.5, font.main=4)
	#title(paste0("iHS ", pop1), adj=0.01, font.main=2)
	#manhattanplot(wg.ihs.pop1, pval= F, pch= 20)
	#title(paste0("iHS score for ",pop1), adj= 0.01, font.main=2)
	#dev.off()
#}


### XPEHH for pop1 vs. pop2
#Populations_List1 = c("BurkinaFaso_Fulani_Banfora","BurkinaFaso_Fulani_Tindangou","BurkinaFaso_Fulani_Ziniare","Cameroon_Fulani_Tcheboua","Chad_Fulani_Bongor","Chad_Fulani_Linia","Guinea_Fulani","Mali_Fulani_Diafarabe","Mali_Fulani_InnerDelta","Mauritania_Fulani_Assaba","Niger_Fulani_Abalak","Niger_Fulani_Ader","Niger_Fulani_Balatungur","Niger_Fulani_Diffa","Niger_Fulani_Zinder","Senegal_Fulani_Linguere")
Populations_List1 = c("Wodaabe","Eastern", "Western")
Populations_List2 = c("Senegal_Bedik","Chad_Toubou","Yemen_Yemeni")


for (p2 in 1:length(Populations_List2)) {
	pop2=Populations_List2[[p2]]
	#wg.res.pop2 <- readRDS(file= paste0("rehh_outputs/",pop2,"_resdata.rds"))
	
	for (p1 in 1:length(Populations_List1)) {
		pop1=Populations_List1[[p1]]
		#wg.res.pop1 <- readRDS(file= paste0("rehh_outputs/",pop1,"_resdata.rds"))
		
		
		#wg.xpehh <- ies2xpehh(wg.res.pop1,wg.res.pop2, p.adjust.method= "BH")
		#saveRDS(wg.xpehh, file= paste0("rehh_outputs_adjusted-Pvalue/",pop1,"_vs_",pop2,"_xpehh_data.rds"))
		wg.xpehh <- readRDS(file= paste0("rehh_outputs_adjusted-Pvalue/",pop1,"_vs_",pop2,"_xpehh_data.rds")) 
		# If you saved the data you don't need to re-calculate wg.xpehh again.
		

		# First option to calculate the threshold
		# thres <- unname(15*(quantile(wg.ihs.pop1$ihs$LOGPVALUE, prob=0.99, na.rm=T)))
		# cat('\n',pop1,'XPEHH',thres,'\n')
		# 
		# Second option to predefine the threshold
		thres <- 5
		cat('\n',pop1,'XPEHH',thres,'\n')
		#
		cr_wg.xpehh <- calc_candidate_regions(wg.xpehh, ignore_sign= T, window_size= 1E4, threshold= thres, pval= T)
		cr_wg.xpehh
		genes <- find_genes(cr)
		paste(genes,collapse=' ')
		# write.table(genes,file=paste0("list_of_genes_",pop1,"_vs_",pop2,".csv"))
		# write.table(wg.xpehh,file=paste0("wg.xpehh_",pop1,"_vs_",pop2,".csv"))

		### XPEHH pop1 vs. pop2
		pdf(file= paste0("rehh_Plots/Edited_XPEHH_plot/xpehh.plot_",pop1,"_vs_",pop2,".pdf"),height=8,width= 14)
		# png(paste0("rehh_Plots/Edited_XPEHH_plot/xpehh_",pop1,"_vs_",pop2,".png"), height=30, width=40, units='cm', res=300)
		par(mfrow=c(2,1))
		manhattanplot(wg.xpehh, pval= T, threshold= c(-thres,thres),sub= paste(genes, collapse=" "), cex= 0.6, adj= 0.5, font.main=4,
		cr=cr_wg.xpehh, cr.lab.cex=0)
		title(paste0("Adjusted p-values of ",pop1," vs ",pop2), adj= 0.01, font.main=2)
		manhattanplot(wg.xpehh, pval= F, adj= 0.5, font.main=4)
		title(paste0("XPEHH score for ",pop1," vs ",pop2), adj= 0.01, font.main=2)
		dev.off()
		
		## Rsb pop1 vs. pop2. Rsb is really opcional because the results will be very similar than XPEHH.
		# pdf(file= paste0("rehh_Plots/Edited_Rsb_plot/Rsb.plot_",pop1,"_vs_",pop2,".pdf"),height= 5,width= 14)
		# png(paste0("rehh_Plots/Edited_Rsb_plot/Rsb_",pop1,"_vs_",pop2,".png"), height= 5, width= 14, units= 'in', res= 300)
		# manhattanplot(wg.rsb, pval= T,  main= paste0("Rsb score for ",pop1," vs ",pop2), threshold= c(-5,5))
		# dev.off()
}
}
