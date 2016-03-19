library(Biostrings)
library(sqldf)

matc <- readDNAStringSet("chr13_m27_maternal.fa")
patc <- readDNAStringSet("chr13_m27_paternal.fa")
myexons <- read.table("c13_exons.txt",header=T)


getCDSStringSet <- function(chr,exons) {
	genelist <- sqldf("SELECT Gene_ID from exons where direction='reverse' group by Gene_ID")[,1]
	ir <- IRanges(c(exons$Start),c(exons$End), names=exons$Gene_ID)
	v <- Views(unlist(chr),ir)
	dna <- DNAStringSet(v)
	dna <- DNAStringSet(sapply(unique(dna@ranges@NAMES),function(x) c(unlist(dna[dna@ranges@NAMES==x]))))
	dna[genelist] <- reverseComplement(dna[genelist]) 
	return(dna)
}

getFirstStop <- function(aa) {
	return(sapply(aa,function(x)  regexpr("\\*",x)[1]))
}
		

mat_dna <- getCDSStringSet(matc,myexons)
pat_dna <- getCDSStringSet(patc,myexons)

mat_aa <- translate(mat_dna,if.fuzzy.codon="solve")
pat_aa <- translate(pat_dna,if.fuzzy.codon="solve")


mat_stop <- getFirstStop(mat_aa) 
pat_stop <- getFirstStop(pat_aa) 

test <- (mat_aa@ranges@width/mat_stop)
nonstop_mat <- test[test<1]
nonsense_mat <- test[test>1]

test <- (pat_aa@ranges@width/pat_stop)
nonstop_pat <- test[test<1]
nonsense_pat <- test[test>1]

nonstop_mutants_mat <- names(nonstop_mat[!names(nonstop_mat)%in%names(nonstop_pat)])
nonstop_mutants_pat <- names(nonstop_pat[!names(nonstop_pat)%in%names(nonstop_mat)])

nonsense_mutants_mat <- names(nonsense_mat[!names(nonsense_mat)%in%names(nonsense_pat)])
nonsense_mutants_pat <- names(nonsense_pat[!names(nonsense_pat)%in%names(nonsense_mat)])

dna_align <- pairwiseAlignment(mat_dna,pat_dna)
aa_align <- pairwiseAlignment(mat_aa,pat_aa)

setdiff.data.frame <- function(A,B) A[ !duplicated( rbind(B,A) )[ -seq_len(nrow(B))] , ]


lev_codon <- function(codon) {
	x <- as.vector(strsplit(as.character(codon),""))
	codon <- x
	l2 <- 	0 + 
		(codon[[1]][2]=="A") + 
		(codon[[1]][2]=="G" & grepl(codon[[1]][1],c("ACM"))& grepl(codon[[1]][3],c("AGR"))) +
		(codon[[1]][2]=="T" & grepl(codon[[1]][1],c("CTY"))& grepl(codon[[1]][3],c("AGR"))) +
		(codon[[1]][1]=="T" & codon[[1]][3]=="A" & grepl(codon[[1]][2],c("AGR"))) + 
		(codon[[1]][1]=="A" & codon[[1]][2]=="G") +
		(codon[[1]][1]=="T" & codon[[1]][2]=="T") +
		(codon[[1]][1]=="T" & codon[[1]][2]=="G" & grepl(codon[[1]][3],c("CTY"))) +
		(codon[[1]][1]=="A" & codon[[1]][2]=="T" & grepl(codon[[1]][3],c("ACTMWYH")))
	l4 <- 	0 + 	
		(codon[[1]][2]=="C" & grepl(codon[[1]][1],c("ACTMWYH"))) + 
		(grepl(codon[[1]][2],c("GTK")) & grepl(codon[[1]][1],c("CGS")))
	l0 <- 3 - (l2 + l4)
	c(l0,l2,l4)	
}	

lev_nuc <- function(codon) {
	x <- as.vector(lev_codon(codon))
	x <- paste(lev_codon(codon),collapse="")
	y <- switch (x,
		"111"=c(2,0,4),
		"201"=c(0,0,4),
		"300"=c(0,0,0),
		"210"=c(0,0,2),
		"120"=c(2,0,2)
		)
	if(codon=="TGA"){y<-c(0,2,0)}
	if(codon=="TAA"){y<-c(0,2,2)}
	return(y)
}


synnonlist <- function(dna_a,aa_a,md,pd) {
	mt1 <- mismatchTable(dna_a)[,c(2,4,8,1)]
	mt2 <- mismatchTable(aa_a)[,c(2,4,8,1)]
	mt1$gene <- md@ranges@NAMES[mt1[,4]]
	mt2$gene <- md@ranges@NAMES[mt2[,4]]	
	testsql <- sqldf("select
			mt1.gene,  
			mt1.PatternStart,
			mt1.PatternSubstring,
			mt1.SubjectSubstring,
			mt2.PatternStart as AAStart,
			mt2.PatternSubstring as AA_Mat,
			mt2.SubjectSubstring as AA_Pat			  				
		        from mt1 inner join mt2 on mt1.gene = mt2.gene
		")
	testsql$codon_start <- ceiling(testsql$PatternStart/3)*3-2
	nonsyn <- testsql[(ceiling(testsql[,2]/3))==testsql[,5],]
	ambig <- nonsyn[(as.logical(duplicated(nonsyn[,c(1,8)],fromLast=F)+ duplicated(nonsyn[,c(1,8)],fromLast=T))),]
	syn <- setdiff.data.frame(mt1[,c(5,1:3)],nonsyn[,c(1:4)])
	syn$codon_start <- ceiling(syn$PatternStart/3)*3-2
	nonsyn <- setdiff.data.frame(nonsyn,ambig)
	nonsyn <- nonsyn[,-5]
	ambig <- unique(ambig[,c(1,8)])
	ambig$type <- "nonsyn"

	syn_ambig <- unique(syn[as.logical(duplicated(syn[,c(1,5)])+(duplicated(syn[,c(1,5)],fromLast=T))),])
	syn <- setdiff.data.frame(syn,syn_ambig)
	
	syn_ambig <- syn_ambig[,c(1,5)]
	syn_ambig$type <- "syn"
	ambig <- rbind(ambig,syn_ambig)

	minifunc <- function(X,Y) {
		ir <- IRanges(c(X$codon_start),c(X$codon_start+2), names=X$gene)
		x <- Y[ir@NAMES]
		return(as.character(DNAStringSet(x,start=ir@start,width=ir@width)))
	}

	ambig$mat <- minifunc(ambig,md)
	ambig$pat <- minifunc(ambig,pd)
	nonsyn$mat <- minifunc(nonsyn,md)
	nonsyn$pat <- minifunc(nonsyn,pd)
	syn$mat <- minifunc(syn,md)
	syn$pat <- minifunc(syn,pd)
	
	nonsyn <- nonsyn[,c(1,7:9,5,6)]
	list(synonymous=syn,nonsynonymous=nonsyn,ambiguous=ambig)
}

mysynnon <- synnonlist(dna_align,aa_align,mat_dna,pat_dna)



#####old#####

#mysynnon$syn$gene <- mat_dna@ranges@NAMES[mysynnon$syn$PatternId]
#mysynnon$nonsyn$gene <- mat_dna@ranges@NAMES[mysynnon$nonsyn$PatternId]

#mysynnon$nonsyn$codon_start <- ceiling(mysynnon$nonsyn[,1:4]$PatternStart/3)*3-2
#mysynnon$nonsyn$codon_end <- ceiling(mysynnon$nonsyn[,1:4]$PatternStart/3)*3

synnonlist <- function(dna_a,aa_a) {
#	mt1 <- mismatchTable(dna_a)[,c(2,4,8,1)]
#	mt2 <- mismatchTable(aa_a)[,c(2,4,8,1)]	
#	testsql <- sqldf("select * from mt1 inner join mt2 on mt1.PatternId = mt2.PatternId ")
#	nonsyn <- testsql[(ceiling(testsql[,1]/3))==testsql[,5],]
#	#testmerge <- merge(mt1,mt2,by.x="PatternId", by.y="PatternId")
#	#nonsyn <- testmerge[(ceiling(testmerge[,2]/3))==testmerge[,5],]
#	nonsyn <- nonsyn[,c(1:3,AAPattern=6,AASubject=7,8)]
#	syn <- setdiff.data.frame(mt1,nonsyn[,1:4])
#	list(syn=syn,nonsyn=nonsyn)
#}

#q_f <- myexons[myexons$direction=="forward",]
#q_r <- myexons[myexons$direction=="reverse",]

#genes <- sqldf("select min(Start) as start, max(End) as end, Gene_ID, direction from myexons group by Gene_ID")
#g_f <- sqldf("select min(Start) as start, max(End) as end, Gene_ID, direction from myexons where direction='forward' group by Gene_ID")
#g_r <- sqldf("select min(Start) as start, max(End) as end, Gene_ID, direction from myexons where direction='reverse' group by Gene_ID")



junk_func <- function(dna1,dna2,aa1,aa2) {
	t1 <- pairwiseAlignment(dna1,dna2)
	t2 <- pairwiseAlignment(aa1,aa2)
	mt1 <- mismatchTable(t1)[,c(2,4,8,1)]
	mt2 <- mismatchTable(t2)[,c(2,4,8,1)]	
	syn <- mt1[(mt1[,4]%in%mt2[,4])&(!ceiling(mt1[,1]/3) %in% mt2[,1]),]
	nonsyn <- mt1[!mt1[,1]%in%syn[,1],]
	list(syn=syn,nonsyn=nonsyn)
}

test <- junk_func(mat_dna[1:100],pat_dna[1:100],mat_aa[1:100],pat_aa[1:100])



g_f <- getFirstStop(mat_aa_f,g_f,"mat_stop")
g_f <- getFirstStop(pat_aa_f,g_f,"pat_stop")
g_r <- getFirstStop(mat_aa_r,g_r,"mat_stop")
g_r <- getFirstStop(pat_aa_r,g_r,"pat_stop")

diff_f2 <- g_f[g_f[,5]-g_f[,6]!=0,]
diff_r2 <- g_r[g_r[,5]-g_r[,6]!=0,]

tf_m <- mat_aa_f[diff_f2$Gene_ID]
tf_p <- pat_aa_f[diff_f2$Gene_ID]
tr_m <- mat_aa_r[diff_r2$Gene_ID]
tr_p <- pat_aa_r[diff_r2$Gene_ID]

align_f <- pairwiseAlignment(tf_m,tf_p)
align_r <- pairwiseAlignment(tr_m,tr_p)

diff_f2 <- cbind(diff_f2,match=nmatch(align_f),mis=nmismatch(align_f),mat=as.character(tf_m),pat=as.character(tf_p))
diff_r2 <- cbind(diff_r2,match=nmatch(align_r),mis=nmismatch(align_r),mat=as.character(tr_m),pat=as.character(tr_p))

nonsense <- rbind(diff_f2,diff_r2)

write.table(nonsense,"c13_m27_nonsense.txt",sep="\t",row.names=F,quote=F)

#test[exons$Gene_ID[exons$direction[unique(exons$Gene_ID)]=="reverse"]

#tt <- exons[unique(exons$Gene_ID),]$Gene_ID[exons[unique(exons$Gene_ID),]$direction=="reverse"]
t1 <- exons[unique(exons$Gene_ID),c(3,5)]
t2 <- t1[t1$direction=="reverse",1]

#tv <- exons[unique(exons$Gene_ID),]$Gene_ID[exons[unique(exons$Gene_ID),]$direction=="forward"]
#	ir@metadata <- as.list(exons$direction)	
#	dna[unlist(dna@ranges@metadata)=="reverse"] <- reverseComplement(dna[unlist(dna@ranges@metadata)=="reverse"]) 
#test[tt] <- reverseComplement(dna[tt]) 
#dna[unlist(dna@ranges@metadata)=="reverse"] <- 
