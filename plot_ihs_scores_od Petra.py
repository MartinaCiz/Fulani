import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="whitegrid")
import pylab
import matplotlib.ticker as ticker
from itertools import cycle, islice
import plotly.express as px
import plotly.graph_objects as go
import plotly
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="whitegrid")
import pylab
import matplotlib.ticker as ticker
from itertools import cycle, islice
from plotly.subplots import make_subplots
import json
import cPickle as pickle
#with open(r"d.pck", "rb") as inp:
#	c = pickle.load(inp)
#print c
#quit()


def parse_path_genes():
	path_genes_dict = {}
	file = open("path_genes.csv").read().splitlines()
	for line in file:
		lline = line.split(",")
		path_genes_dict.update({lline[0]:lline[1].split()})
	return path_genes_dict


def todict(string):
	dic ={}
	itms = string.split(";")
	for it in itms:
		it = it.split("=")
		dic.update({it[0]:it[1]})
	
	if "gene" in dic:
		return dic["gene"]
	else:
		return "intergenic"

def get_chrom(transcript):
	chro = transcript.split("_")[1].split(".")[0].lstrip("0")
	return chro 


#refseq = pd.read_csv("GRCh37_latest_genomic.gff", sep="\t",comment ="#")
#refseq["dict"]= refseq["ANNOT"].apply(todict)

#refseq.to_json("refseq.json")

#pops = ["Daza","Maba","Daju","Dangaleat","Kawalib","Zaghawa_Sudan","Bedik","Arabs_Rashayda"]
pops = [
"Fulani_Fouta_Djallon",
"Fulani_Halpularen",
"Bedik",
"Dangaleat",
#"Arabs_Chad",
"Arabs_Baggara",
"Daza",
"Daju",
"Maba",
#"Kabbabish/Kawahla",
"Kawalib",
"Arabs_Rashayda",
"Zaghawa_Sudan"
]

with open("refseq.pck") as input:
	refseq = pickle.load(input)

refseq = refseq[refseq["TYPE"]=="gene"]
refseq["CHROM"]= refseq["TRANSCRIPT"].apply(get_chrom)

path_genes_dict = parse_path_genes()

#refseq = pd.read_json("refseq.json")
#print refseq
#quit()

def gene_in_path(genes, path, path_genes_dict):
	#col = "#D0D0D0"
	genes = genes.split(",")
	inpath = "OUT"
	for gene in genes:
		if gene in path_genes_dict["Immune system"]:
			#col = "#E31F1F"
			inpath = "IMMUNE"
		if gene in path_genes_dict["Immune disease"]:
			#col = "#E31F1F"
			inpath = "IMMUNE"
		if gene in path_genes_dict["Carbohydrate metabolism"]:
			#col = "#E31F1F"
			inpath = "CARBOHYDRATE"
		if gene in path_genes_dict["Lipid metabolism"]:
			#col = "#E31F1F"
			inpath = "LIPID"
		if gene in path_genes_dict["Infectious disease: bacterial"]:
			#col = "#E31F1F"
			inpath = "INFECT"
		if gene in path_genes_dict["Infectious disease: viral"]:
			#col = "#E31F1F"
			inpath = "INFECT"
		if gene in path_genes_dict["Infectious disease: parasitic"]:
			#col = "#E31F1F"
			inpath = "INFECT"			
		
		if gene in path_genes_dict["Amino acid metabolism"]:
			#col = "#E31F1F"
			inpath = "AMINO"
		if gene in path_genes_dict["Digestive system"]:
			#col = "#E31F1F"
			inpath = "DIGEST"	
		
	return inpath


def get_annot(chrompos):
	global refseq
	chrom = int(chrompos.split("_")[0])
	pos = int(chrompos.split("_")[1])
	#print chrom, pos,
	ref = refseq[refseq["CHROM"]==str(chrom)]
	#print ref
	df = ref[(ref["START"]<=pos) & (ref["STOP"]>=pos)]
	annot = df["dict"]
	#print ",".join(annot)
	return ",".join(annot)

chromosomes = range(1,23)

#df_in = pd.read_csv(sys.argv[1], sep=",")#, names=['chrom','id','pos','cosi','fst'])
#df_in = df_in.sort_values(["FST"], ascending=False)


#df_in =df_in[df_in.FST>0.05]
#df_in = df_in.head(10000)
#print df_in["Consequence"]
#quit()
#df_in=df_in[df_in['Consequence'].str.contains('|'.join(["frame","stop"]))]
fig = make_subplots(rows=len(pops), cols=22, shared_yaxes=True,horizontal_spacing =0,vertical_spacing=0.05, row_titles = pops) #subplot_titles = chromosomes)

#print df_in[df_in['FST']>0.05]
#quit()


row_cntr = 0
for pop in pops:
	row_cntr += 1
	odd=-1
	cntr=1
	for chrom in range(1,23):
		odd = odd*-1
		if odd==1:
			color = '#A52A2A'
		else:
			color = '#33336B'	
		df=pd.read_csv(pop+"_sahel_tgen_QC_chr"+str(chrom)+".selhaps.hap.gz.ihs.ihs.out.100bins.norm.gz",sep="\t", names=["id","physpos","freq","ihh1","ihh0","ihs","ihs_norm","?"])
		df["CHROMPOS"]= str(chrom)+"_"+df.physpos.astype(str)
		df["ihs_abs"]=df["ihs_norm"].abs()
		#print df.ihs_abs
		df = df[df["ihs_abs"]>3.5]
		df["GENE"]= df["CHROMPOS"].apply(get_annot)
		df["INPATH"] = df["GENE"].apply(gene_in_path, args=("Immune system", path_genes_dict))
		df_immune = df[df["INPATH"]=="IMMUNE"]
		df_carbo = df[df["INPATH"]=="CARBOHYDRATE"]
		df_lipid = df[df["INPATH"]=="LIPID"]
		df_infect = df[df["INPATH"]=="INFECT"]
		df_amino = df[df["INPATH"]=="AMINO"]
		df_digest = df[df["INPATH"]=="DIGEST"]
		df_out = df[df["INPATH"]=="OUT"]
		#print "IN"
		#print df_in["GENE"]
		#print "OUT"
		#print df_out["GENE"]
		#quit()
#		fig.append_trace(go.Scatter(hovertext = df.GENE, name=chrom, x=df.physpos.astype(int), y=df.ihs_abs.astype(float),marker=dict(color=color),mode='markers'), col = cntr, row=row_cntr)
		fig.append_trace(go.Scatter(hovertext = df_out.GENE, name=chrom, x=df_out.physpos.astype(int), y=df_out.ihs_abs.astype(float),marker=dict(color= "#B6AEAE"),mode='markers'), col = cntr, row=row_cntr) #gray

		fig.append_trace(go.Scatter(hovertext = df_immune.GENE, name=chrom, x=df_immune.physpos.astype(int), y=df_immune.ihs_abs.astype(float),marker=dict(color= "#DC0808"),mode='markers'), col = cntr, row=row_cntr) #red
		fig.append_trace(go.Scatter(hovertext = df_carbo.GENE, name=chrom, x=df_carbo.physpos.astype(int), y=df_carbo.ihs_abs.astype(float),marker=dict(color= "#0000FF"),mode='markers'), col = cntr, row=row_cntr) #navy
		fig.append_trace(go.Scatter(hovertext = df_lipid.GENE, name=chrom, x=df_lipid.physpos.astype(int), y=df_lipid.ihs_abs.astype(float),marker=dict(color= "#FFFF00"),mode='markers'), col = cntr, row=row_cntr) #yellow
		fig.append_trace(go.Scatter(hovertext = df_infect.GENE, name=chrom, x=df_infect.physpos.astype(int), y=df_infect.ihs_abs.astype(float),marker=dict(color= "#EB4B67"),mode='markers'), col = cntr, row=row_cntr) #pink
		fig.append_trace(go.Scatter(hovertext = df_amino.GENE, name=chrom, x=df_amino.physpos.astype(int), y=df_amino.ihs_abs.astype(float),marker=dict(color= "#1CB530"),mode='markers'), col = cntr, row=row_cntr) #light green
		fig.append_trace(go.Scatter(hovertext = df_digest.GENE, name=chrom, x=df_digest.physpos.astype(int), y=df_digest.ihs_abs.astype(float),marker=dict(color= "#00ADE6"),mode='markers'), col = cntr, row=row_cntr) #light blue

		
		cntr += 1

for row in range(1,len(pops)+1):
	odd = -1
	for chr in range(1,23):
		odd = odd*-1
		if odd==1:
			color = '#A52A2A'
		else:
			color = '#33336B'	
		fig.update_xaxes(title_text= str(chr), row=row, col=chr, color = color)

fig.update_layout(showlegend=False, xaxis=dict(showgrid=False),yaxis=dict(showgrid=False), height = 2500, title = {"yanchor":"bottom"})
fig.update_xaxes(showgrid=False, zeroline=False, showticklabels=False)
fig.update_yaxes(showgrid=False, zeroline=False)


#fst.plot.scatter(x="POS", y="FST", c="b")

fig.show()
#fig.write_image("fst".jpeg",width=1800, height=900,scale=2)
