#load("regulonDB_list.RData")
#regdb.net <- regulonDB_list

regdb.net=read.delim(gzfile("network_tf_gene_ecoli_regdb_isstrong.tsv.gz"),head=F)
regdb.regulons=tapply(as.character(regdb.net$V5),as.character(regdb.net$V3),function(i)unique(i))

## Try it after grouping all genes based on their unique combination of regulators
tmp <- tapply(as.character(regdb.net$V3),as.character(regdb.net$V5),function(i)paste(sort(unique(i)),collapse=' '))
tmp <- data.frame(gene=names(tmp), tfs=tmp)
regdb.combi.regulons=tapply(as.character(tmp$gene),as.character(tmp$tfs),function(i)unique(i))

if ( ! exists( 'use.combi' ) ) use.combi <- TRUE
if ( use.combi ) { ## combi-regs
    regdb.regulons <- regdb.combi.regulons
}

regdb.regulons <- regdb.regulons[ which( sapply( regdb.regulons, length ) > 2 ) ]

all.regdb.genes <- unique( unlist( regdb.regulons ) )
all.cm.genes <- unique( unlist( lapply( e$clusterStack, '[[', 'rows' ) ) )

require( parallel )
pvals <- mclapply( 1:e$k.clust, function( k ) {
    #print(k)
    rows <- tolower(e$clusterStack[[ k ]]$rows)
    out <- do.call( rbind, lapply( names(regdb.regulons), function( tf ) {
        reg <- regdb.regulons[[ tf ]]
        if ( length( rows ) <= 2 || length( reg ) < 2 || sum(rows%in%reg) < 2 ) return( NULL )
        out <- phyper( sum(rows%in%reg), length(reg), length(all.regdb.genes)-length(reg), length(rows), lower=F, log=T )
        out <- data.frame( k=k, tf=tf, nrow=length(rows), nreg=length(reg), n=sum(rows%in%reg), pv=out )
    } ) )
} )

pvals <- do.call( rbind, pvals )
pvals$fdr <- log( p.adjust( exp(pvals$pv), 'BH' ) )
pvals$qv <- log( p.adjust( exp(pvals$pv), 'bonferroni' ) )

#tmp <- subset( pvals, n > 2 & qv <= log(0.01) )
#print( length( unique( tmp$k ) ) )
#print( length( unique( tmp$tf ) ) )

tmp <- subset( pvals, n > 2 & fdr <= log(0.01) )
#print( length( unique( tmp$k ) ) )
#print( length( unique( tmp$tf ) ) )

eco.reg.hits <- tmp

## try rand index -- see https://en.wikipedia.org/wiki/Rand_index

## pairsA <- unique( unlist( lapply( e$clusterStack, function(i) if(length(i$rows)>1)
##                                  apply(t(combn(i$rows,2)),1,paste,collapse=' ') else NULL ) ) )
## pairsB <- unlist( lapply( regdb.regulons[sapply(regdb.regulons,length)>=2],
##                          function(i) apply(t(combn(i,2)),1,paste,collapse=' ') ) )
## allPairs <- apply(t(combn(unique(unlist( regdb.regulons[sapply(regdb.regulons,length)>=2] )), 2)),
##                   1,paste,collapse=' ')
## R <- ( sum(pairsA %in% pairsB) + sum(!allPairs%in%c(pairsA, pairsB)) ) / length(allPairs)

require( flexclust )

tab <- sapply( e$clusterStack, function(i)
              sapply( regdb.regulons[sapply(regdb.regulons,length)>=2], function(j)
                     sum(unique(tolower(i$rows)) %in% j) ) )
R <- randIndex( as.table(tab), correct=T )
R1 <- randIndex( as.table(tab), correct=F )

