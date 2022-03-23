# First create folders in your working directory
# mkdir rehh_outputs rehh_Tables rehh_Plots rehh_Plots/QQ_plot rehh_Plots/XPEHH_plot rehh_Plots/Rsb_plot

# rehh script to run in R

library(rehh)
library(qqman)

# Full list
# Populations_List = c("BurkinaFaso_Fulani_Banfora","BurkinaFaso_Fulani_Tindangou","BurkinaFaso_Fulani_Ziniare","BurkinaFaso_Gurmantche","BurkinaFaso_Gurunsi","BurkinaFaso_Mossi","Cameroon_Fulani_Tcheboua","Chad_Daza","Chad_Fulani_Bongor","Chad_Fulani_Linia","Chad_Kanembu","Chad_Laal","Chad_Sara","Chad_Toubou","Guinea_Fulani","Lebanon_LebaneseChristian","Lebanon_LebaneseDruze","Lebanon_LebaneseMuslim","Mali_Fulani_Diafarabe","Mali_Fulani_InnerDelta","Mauritania_Fulani_Assaba","Niger_Fulani_Abalak","Niger_Fulani_Ader","Niger_Fulani_Balatungur","Niger_Fulani_Diffa","Niger_Fulani_Zinder","Senegal_Bedik","Senegal_Fulani_Linguere","Senegal_Halpularen","Sudan_Daju","Sudan_NubaKoalib","Sudan_Nubian","Sudan_Zaghawa","Yemen_Yemeni")

#pooled pops
Populations_List = c("Wodaabe", "Eastern", "Western")


for (p in 1:length(Populations_List)) {
	pop1=Populations_List[[p]]
	
	for(i in 1:22){
	hap_file= paste("Populations_inputs/",pop1,"/",pop1,"_",i,"_hapguess_switch.out",sep= "")
	mapFile= paste("Genetic_map/this_DB_chr",i,"_map.inp",sep= "")
	data_pop1 <- data2haplohh(hap_file= hap_file, map_file= mapFile, chr.name= i, allele_coding= "map")
	res.pop1 <- scan_hh(data_pop1, maxgap= 200000)
	if(i==1){wg.res.pop1 <-res.pop1}else{wg.res.pop1<-rbind(wg.res.pop1,res.pop1)}
	}
	
	# wg.ihs.pop1 <- ihh2ihs(wg.res.pop1, freqbin= 0.05, min_maf= 0.05)
	wg.ihs.pop1 <- ihh2ihs(wg.res.pop1, freqbin= 0.05, min_maf= 0.05, p.adjust.method= "BH") # P-values are adjusted using the Benjamini & Hochberg (1995) correction

	saveRDS(wg.res.pop1, file= paste0("rehh_outputs/",pop1,"_resdata.rds"))
	saveRDS(wg.ihs.pop1, file= paste0("rehh_outputs/",pop1,"_ihsdata.rds"))
	
	ihs_pop1 <- na.omit(wg.ihs.pop1$ihs)
	map <- data.frame(ID=mapFile[1], POSITION=mapFile[3], Anc=mapFile[4], Der=mapFile[5])
	ihsMerge <- merge(map, ihs_pop1, by= "POSITION")
	# signals <- ihsMerge[ihsMerge[,7]>=0.05]
	# sigpos <- signals[,4]
	
	write.table(ihs_pop1, file= paste0("rehh_Tables/",pop1,"_iHSresult.txt"), col.names=T, row.names=T, quote=F, sep="\t")
	# write.table(signals, file= paste0("rehh_Tables/",pop1,"_sigOut.txt"), col.names=T, row.names=F, quote=F, sep="\t")
	write.table(wg.ihs.pop1$frequency.class, file= paste0("rehh_Tables/",pop1,"_iHSfrq.txt"), col.names=T, row.names=F, quote=F, sep="\t")
	
	# wg.ihs.pop1 <- readRDS(file= paste0("rehh_outputs/",pop1,"_ihsdata.rds"))
	# cr <- calc_candidate_regions(wg.ihs.pop1, ignore_sign= T, window_size= 1, threshold= 7)
	
	# pdf(file= paste0("rehh_Plots/",pop1,"_iHS.plot.pdf"), height= 4, width= 8)
	png(paste0("rehh_Plots/",pop1,"_iHS.plot.png"), height= 5, width= 14, units= 'in', res= 300)
	manhattanplot(wg.ihs.pop1, main= paste0("iHS ", pop1), pval= T, pch= 20, threshold= c(-5,5)) # , threshold = c(-thresh, thresh),  cr= cr_pop1, threshold= c(-1.5,1.5), ylim= c(-2.5,2.5) 
	dev.off()
	
	# Gaussian Distribution and Q-Q plots
	IHS <- wg.ihs.pop1$ihs[["IHS"]]
	png(paste0("rehh_Plots/QQ_plot/QQ_",pop1,".png"), height= 700, width= 440, res= NA, units= "px", type= "cairo")
	layout(matrix(1:2,2,1))
	distribplot(IHS, main="iHS", qqplot = F)
	distribplot(IHS, main="iHS", qqplot = T)
	dev.off()
	
}
"Done"

# Populations_List1 = c("BurkinaFaso_Fulani_Banfora","BurkinaFaso_Fulani_Tindangou","BurkinaFaso_Fulani_Ziniare","Cameroon_Fulani_Tcheboua","Chad_Fulani_Bongor","Chad_Fulani_Linia","Guinea_Fulani","Mali_Fulani_Diafarabe","Mali_Fulani_InnerDelta","Mauritania_Fulani_Assaba","Niger_Fulani_Abalak","Niger_Fulani_Ader","Niger_Fulani_Balatungur","Niger_Fulani_Diffa","Niger_Fulani_Zinder","Senegal_Fulani_Linguere")
# Populations_List2 = c("Senegal_Bedik","Chad_Toubou","Yemen_Yemeni")

# for (p2 in 1:length(Populations_List2)) {
	# pop2=Populations_List[[p2]]
	# wg.res.pop2 <- readRDS(file= paste0("rehh_outputs/",pop2,"_resdata.rds"))
	
	# for (p1 in 1:length(Populations_List1)) {
		# pop1=Populations_List[[p1]]
		
		# wg.res.pop1 <- readRDS(file= paste0("rehh_outputs/",pop1,"_resdata.rds"))	
		
		### XPEHH pop1 vs. pop2
		# wg.xpehh<-ies2xpehh(wg.res.pop1,wg.res.pop2)
		# pdf(file= paste0("rehh_Plots/XPEHH_plot/xpehh.plot_",pop1,"_vs_",pop2,".pdf"),height= 5,width= 14)
		# png(paste0("rehh_Plots/XPEHH_plot/xpehh_",pop1,"_vs_",pop2,".png"), height= 5, width= 14, units= 'in', res= 300)
		# manhattanplot(wg.xpehh, pval= T, main= paste0("XPEHH ",pop1," vs ",pop2), threshold= c(-5,5))
		# dev.off()
		
		### Rsb pop1 vs. pop2
		# wg.rsb <- ines2rsb(wg.res.pop1,wg.res.pop2)
		# pdf(file= paste0("rehh_Plots/Rsb_plot/Rsb.plot_",pop1,"_vs_",pop2,".pdf"),height= 5,width= 14)
		# png(paste0("rehh_Plots/Rsb_plot/Rsb_",pop1,"_vs_",pop2,".png"), height= 5, width= 14, units= 'in', res= 300)
		# manhattanplot(wg.rsb, pval= T,  main= paste0("Rsb ",pop1," vs ",pop2), threshold= c(-5,5))
		# dev.off()
# 	}
# }
