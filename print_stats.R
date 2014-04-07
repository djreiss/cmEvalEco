print(sprintf('k =         %d', e$k.clust))
print(sprintf('Avg n.row = %.2f', mean(sapply(e$clusterStack,function(cl)length(cl$rows)))))
print(sprintf('Avg resid = %.2f', mean(sapply(e$clusterStack,function(cl)cl$resid))))
print(sprintf('Avg motifs with E < 1e-3 per cluster = %.2f', mean(sapply(e$clusterStack,function(cl)
                                                          sum(cl$e.val <= 1e-3, na.rm=T)))))
print(sprintf('Avg motifs with E < 1 per cluster = %.2f', mean(sapply(e$clusterStack,function(cl)
                                                          sum(cl$e.val <= 1, na.rm=T)))))
print(sprintf('Avg motifs with E < 10 per cluster = %.2f', mean(sapply(e$clusterStack,function(cl)
                                                          sum(cl$e.val <= 10, na.rm=T)))))

print(sprintf('Fraction of clusters with a E < 1e-3 motif = %.2f',
              mean(sapply(e$clusterStack,function(i)min(i$e.val))<=1e-3,na.rm=T)))
print(sprintf('Fraction of clusters with a E < 1 motif = %.2f',
              mean(sapply(e$clusterStack,function(i)min(i$e.val))<=1,na.rm=T)))
print(sprintf('Fraction of clusters with a E < 10 motif = %.2f',
              mean(sapply(e$clusterStack,function(i)min(i$e.val))<=10,na.rm=T)))

c.per.g <- rep( 0, max(4203,nrow(e$ratios$ratios)) )
names( c.per.g )[ 1:nrow(e$ratios$ratios) ] <- rownames(e$ratios$ratios)
tab <- table(unlist(lapply(e$clusterStack,'[[', 'rows')))
c.per.g[ names(tab) ] <- tab
print(sprintf('Avg clusters per gene = %.2f', mean(tab)))
print(sprintf('Fraction of all genes included = %.2f', mean(c.per.g>0)))

if ( exists( 'eco.hits2' ) ) {
  for ( max.motifs in c( 402, Inf ) ) {
    ##max.motifs <- Inf ##402 ## Maximum # of motifs per clustering to allow (same as # of motif clusters in eco ensemble)
    eval.cutoff <- Inf
    if ( ! is.infinite( max.motifs ) ) {
      if ( sum( ! is.na(e.vals) ) > max.motifs ) eval.cutoff <- sort( e.vals[!is.na(e.vals)] )[ max.motifs ]
    }
    tmp <- subset( eco.hits2, distance <= 0.99 ) ## we used 0.99 for the cutoff for the ensemble analysis, but that also had a p-value and frac cutoff so use 0.98 here
    tmp <- tmp[ order( tmp$distance ), ]
    b.m <- t( apply( sapply( strsplit( as.character( tmp$motif ), '_' ), '[', 2:3 ), 2, as.integer ) )
    tmp <- tmp[ e.vals[ b.m ] <= eval.cutoff, ] # use 100 as e-value cutoff (we did not use this for ensemble analysis; what if we did?)
    tmp <- subset( tmp, ! duplicated( tmp$motif ) )
    print( sprintf('Total unique RegulonDB motifs (out of %d) detected (max motifs %.0f) = %d',
                   length(m.eco2), max.motifs, length( unique( tmp$eco.tf ) )) )

    clusts <- as.integer( sapply( strsplit( as.character( tmp$motif ), '_' ), '[', 2 ) )
    n.hits <- rep( 0, e$k.clust ); names( n.hits ) <- 1:e$k.clust
    clusts.tab <- table( clusts )
    n.hits[ names( clusts.tab ) ] <- clusts.tab
    print( sprintf('Avg number of RegulonDB motifs detected per cluster (max motifs %.0f) = %.2f',
                   max.motifs, mean( n.hits ) ) )
  }
}


if ( exists( 'eco.reg.hits' ) ) {
    print( sprintf('Total RegulonDB regulons detected (max regulons %.0f) = %d',
                   sum(sapply(regdb.regulons,length)>=2), length( unique( eco.reg.hits$tf ) )) )
    print( sprintf('Total clusters matching a regulon = %d',length( unique( eco.reg.hits$k ) )) )
}

if ( exists( 'R1' ) ) {
    print( sprintf('RAND index vs. regulons (corrected)   = %.4f', R ) )
    print( sprintf('RAND index vs. regulons (uncorrected) = %.4f', R1 ) )
}
