########## Try to do the same to match up motifs w/ sigma factors -- use motifDBs/RegulonDB/RegulonDATADist/promoter.txt
tmp <- readLines( 'motifDBs/RegulonDB/RegulonDATADist/promoter.txt' )
skip <- grep( 'Table name: PROMOTER', tmp, fixed=T ) + 1
rm( tmp )
tab2 <- read.delim( 'motifDBs/RegulonDB/RegulonDATADist/promoter.txt', skip=skip, head=F )
tab2 <- subset( tab2, ! is.na( V4 ) & V5 != '' )
combos <- grep( ',', tab2$V5, fixed=T )
for ( i in combos ) {
  ii <- strsplit( as.character( tab2$V5[ i ] ), ',', fixed=T )[[ 1 ]]
  tab2$V5[ i ] <- ii[ 1 ]
  tab2 <- rbind( tab2, tab2[ i, ] )
  tab2$V5[ nrow( tab2 ) ] <- ii[ 2 ]
}
n.tot3 <- sort( table( as.character( tab2$V5 ) ), decreasing=T )
n.tot3 <- n.tot3[ n.tot3 >= 3 ]
tab2 <- subset( tab2, V5 %in% names( n.tot3 ) )

out$e$dlf( 'PromoterSet.txt', 'http://regulondb.ccg.unam.mx/data/PromoterSet.txt' )
tmp <- readLines( 'PromoterSet.txt' )
skip <- grep( 'Evidence that supports the existence of the promoter', tmp, fixed=T ) + 1
rm( tmp )
tab2 <- read.delim( 'PromoterSet.txt', skip=skip, head=F )
tab2 <- subset( tab2, ! is.na( V4 ) & V5 != '' & V5 != 'unknown' )
## combos <- grep( ',', tab2$V5, fixed=T )
## for ( i in combos ) {
##   ii <- strsplit( as.character( tab2$V5[ i ] ), ',', fixed=T )[[ 1 ]]
##   tab2$V5[ i ] <- ii[ 1 ]
##   for ( j in 2:length( ii ) ) {
##     tab2 <- rbind( tab2, tab2[ i, ] )
##     tab2$V5[ nrow( tab2 ) ] <- ii[ j ]
##   }
## }
n.tot3 <- sort( table( as.character( tab2$V5 ) ), decreasing=T )
n.tot3 <- n.tot3[ n.tot3 >= 3 ]
tab2 <- subset( tab2, V5 %in% names( n.tot3 ) & V4 != 0 )
tab2 <- subset( tab2, V7 != '[ICWHO|W|Inferred computationally without human oversight]' )
tab2 <- subset( tab2, V7 != '[HIPP|W|Human inference of promoter position]' )
tab2 <- subset( tab2, V7 != '[AIPP|W|Automated inference of promoter position]' )

if ( ! exists( 'm.eco3' ) ) {
  m.eco3 <- mclapply( names( n.tot3 ), function( sf ) {
    m.tmp <- lapply( nc, function( i ) rep( FALSE, length=i ) ); names( m.tmp ) <- names( nc )
    scans <- subset( tab2, V5 == sf )
    if ( nrow(scans) <= 0 ) return( lapply( m.tmp, function( i ) integer(0) ) ) ##Matrix( m.tmp ) )
    for ( iii in 1:nrow( scans ) ) {
      if ( scans$V3 == 'forward' ) inds <- scans$V4[ iii ] + (-45):5
      else if ( scans$V3 == 'reverse' ) inds <- scans$V4[ iii ] + (-5):45
      chr <- names( e$genome.info$genome.seqs )[ 1 ] ##levels( scans$Seq )[ scans$Seq[ iii ] ] ## gene
      m.tmp[[ chr ]][ inds ] <- TRUE ##scans$posns[j]:(scans$posns[j]+wi-1),1] = 1
    }
    cat( sf, which( names( n.tot3 ) == sf ), length( n.tot3 ), nrow( scans ), "\n" ) ##sum( m.tmp ), "\n" )
    lapply( m.tmp, which ) ##which( m.tmp ) ##Matrix( m.tmp )
  } )
  names( m.eco3 ) <- names( n.tot3 )
}

if ( ! exists( 'eco.hits3' ) ) {
  tmp <- mclapply( names( m.eco3 ), function( i ) {
    print( i )
    x <- m.eco3[[ i ]][[ 1 ]]
    if ( length( x ) <= 0 ) return()
    tout <- rep( 1, length( m ) )
    for ( j in 1:length( m ) ) {
      y <- m[[ j ]] ## Assumes it was un-listed in new_compare_regulondb.R
      if ( length( y ) <= 0 ) next
      tmp <- x %in% y
      if ( ! any( tmp ) ) next
      tmp <- ( sum( ! tmp ) + sum( ! ( y %in% x ) ) ) / length( unique( c( x, y ) ) ) ##sum( x != y ) / sum( x | y )
      tout[ j ] <- tmp
    }
    tmp <- which( tout < 1 )
    if ( length( tmp ) <= 0 ) return()
    ##write.table( cbind( i, tmp, out[tmp] ), quote=F, sep='\t', row.names=F, col.names=F, file=bzfile( fname ) )
    data.frame( i, names(m)[tmp], tout[tmp] )
  }, mc.preschedule=F )
  eco.hits3 <- do.call( rbind, tmp )
  rm( tmp )
  colnames( eco.hits3 ) <- c( 'eco.sf', 'motif', 'distance' )
}

ttmp4 <- list()
n.tot <- length( unique( unlist( out$motif.clusts[ 1:out$mc.length ] ) ) )
n.tot4 <- sort( table( eco.hits3$eco.sf ), decreasing=T )
for ( i in 1:out$mc.length ) {
  print(i)
  tmp <- subset( eco.hits3, motif %in% out$motif.clusts[[ i ]] & distance <= 0.999 )
  if ( nrow( tmp ) <= 0 ) next
  n <- sort( table( tmp$eco.sf ), decreasing=T )
  nt <- n.tot4[ names( n ) ]
  ttmp4[[ i ]] <- data.frame( p.val=phyper( n, nt - n, n.tot, length(out$motif.clusts[[i]]), lower=F ),
                             frac=n / length( out$motif.clusts[[ i ]] ),
                             n.mots.in.clust=length(out$motif.clusts[[i]]), n.mots.hit=nt )
}
p.vals <- do.call( c, lapply( ttmp4, function( i ) i[,1] ) )
q.vals <- p.adjust( p.vals, method='BH' )
ind <- 1; for ( i in 1:length( ttmp4 ) ) { if ( is.null( ttmp4[[i]] ) || nrow( ttmp4[[i]] ) <= 0 ) next;
                                               ttmp4[[ i ]]$q.val <- q.vals[ ind:(ind+nrow(ttmp4[[i]])-1) ];
                                               ind <- ind + nrow(ttmp4[[i]]) }
for ( i in 1:length( ttmp4 ) ) { if ( is.null( ttmp4[[ i ]] ) ) next; ttmp4[[ i ]] <- ttmp4[[ i ]][ order( ttmp4[[ i ]]$q.val ), ] }
##hits4 <- sapply( ttmp4, function( i ) which( i[,1] * out$mc.length * length( n.tot4 ) <= 0.01 & i[,2] >= 0.1 ) )
hits4 <- sapply( ttmp4, function( i ) which( i$q.val <= 0.001 & i[,2] >= 0.1 ) )
print(length(sort(table(names(unlist(sapply(hits4,function(i)if(length(i)>1)i[1] else i)))))))
print(sort(table(names(unlist(sapply(hits4,function(i)if(length(i)>1)i[1] else i))))))
## Sigma32 Sigma24 Sigma54 Sigma28 Sigma38 Sigma70
##       3       4       5       9      16      32
