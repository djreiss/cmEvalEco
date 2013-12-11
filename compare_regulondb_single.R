require( multicore )
e$parallel.cores <- e$parallel.cores.motif <- 4
options( cores=4 )
debug.on()

p.cutoff <- 1e-6 ##1e-5

## Get the equivalent of 'm' output from 'new.cluster.motifs()' in cmonkey-ensemble-funcs.R
seq.type <- 1 ##'upstream meme'
ms <- e$meme.scores[[ seq.type ]]
widths1 <- do.call( rbind, mclapply( 1:e$k.clust, function( i ) {
  out <- c( 0, 0, 0 )
  if ( is.null( ms[[ i ]] ) || ms[[ i ]] == '' ) return( out )
  mm <- ms[[ i ]]$meme.out
  for ( j in 1:length( mm ) ) if ( ! is.null( mm[[ j ]]$width ) ) out[ j ] <- mm[[ j ]]$width
  out
} ) )
if ( all( widths1[ ,3 ] == 0 ) ) widths1 <- widths1[ ,1:2 ]
motif.widths <- widths1; rm( widths1, ms )

motifs <- unlist( lapply( 1:nrow( motif.widths ), function( i ) {
  ii <- which( motif.widths[ i, ] > 0 )
  if ( length( ii ) <= 0 ) return( NULL )
  paste( 'MOT', i, ii, sep='_' ) } ) )

nc <- sapply( e$genome.info$genome.seqs, nchar )
nc.cumsum <- c( 0, cumsum( nc ) )[ 1:length( nc ) ]; names( nc.cumsum ) <- names( nc )

#if ( ! exists( 'm' ) ) load( 'filehash/new_motif_shadows_1e-06.RData' )
#m <- m[ motifs ]
## Get the equiv. of 'fimo.out' (from the ensemble stuff) for this bicluster set
if ( ! exists( 'fimo.out' ) ) {
sys.source( "cmonkey-motif-other.R", envir=e )
system( "rm -rvf fimo_out" )
dir.create( 'fimo_out' )
seqs.file <- e$my.tempfile( 'fimo_seqs' )
writeLines( paste( paste( ">", names( e$genome.info$genome.seqs ), sep="" ), e$genome.info$genome.seqs, sep="\n" ),
           con=seqs.file )
inds <- c( seq( 1, e$k.clust, by=100 ), e$k.clust )
files <- mclapply( 1:( length( inds ) - 1 ), function( i ) {
  mots.file <- e$all.motifs.to.mast.file( ks=inds[i]:inds[i+1], seq.type=names(e$mot.weights)[1],
                                         e.value.cutoff=Inf, resid.cutoff=Inf )
  out.file <- e$my.tempfile( 'fimo_out', tmpdir='./fimo_out')
  
  ## Assume using customized version of meme_4.3.0 with MAX_MOTIFS in motif.h changed to 24000
  cmd <- sprintf( './progs/fimo --max-stored-scores 9999999999 --text --verbosity 2 %s %s |bzip2 -c >%s.bz2',
                 mots.file, seqs.file, out.file )
  print( cmd )
  out <- system( cmd, intern=T )
  out.file ## save it to tempfiles because storing it all up on each processor takes too much RAM
}, mc.preschedule=F )
unlink( seqs.file )

fimo.out <- system( sprintf( 'bunzip2 -c fimo_out/fimo_out_*.bz2 | awk \'($6<=%s){print}\'',
                            as.character( p.cutoff ) ), intern=T )
fimo.out <- as.data.frame( do.call( rbind, lapply( fimo.out, function( i ) strsplit( i, '\t' )[[ 1 ]] ) ) )
colnames( fimo.out ) <- strsplit( system( sprintf( 'bunzip2 -c %s.bz2 | head -1', files[1] ), intern=T ), '\t' )[[ 1 ]]
system( "rm -rvf fimo_out" )
fimo.out$Start <- as.integer( as.character( fimo.out$Start ) )
fimo.out$Stop <- as.integer( as.character( fimo.out$Stop ) )
fimo.out$`Log-odds` <- as.numeric( as.character( fimo.out$`Log-odds` ) )
fimo.out$`p-value` <- as.numeric( as.character( fimo.out$`p-value` ) )
tmp <- do.call( rbind, strsplit( as.character( fimo.out$Motif ), '_' ) )
fimo.out$bic <- as.integer( tmp[ ,2 ] )
fimo.out$mot <- as.integer( tmp[ ,3 ] )
fimo.out$Strand <- substr( tmp[ ,1 ], 1, 1 )
rm( tmp )
fimo.out$Motif <- NULL
require( data.table )
fimo.out <- as.data.table( fimo.out )
setkey( fimo.out, bic, mot, Seq, Start )
}  ## if not exists('fimo.out')

if ( ! exists( 'm' ) ) {
m.tmp <- lapply( nc, function( i ) rep( FALSE, length=i ) ); names( m.tmp ) <- names( nc )
p.scans <- fimo.out
m <- mclapply( 1:length( motifs ), function( i ) {
  bi <- as.integer( strsplit( motifs[i], '_' )[[ 1 ]][ 2 ])
  mo <- as.integer( strsplit( motifs[i], '_' )[[ 1 ]][ 3 ])
  scans <- p.scans[ J( bi, mo ) ] ##c( bi, bi ), c( mo, -mo ) ) ]
  if ( nrow(scans) <= 0 ) return( lapply( m.tmp, function( i ) integer(0) ) ) ##Matrix( m.tmp ) )
  scans <- scans[ ! is.na( scans$Start ), ] ##posns ), ]
  if ( nrow(scans) <= 0 ) return( lapply( m.tmp, function( i ) integer(0) ) ) ##Matrix( m.tmp ) )
  for ( iii in 1:nrow( scans ) ) { ##posn in scans$posns ) { ##1:nrow( scans ) ) {
    inds <- scans$Start[ iii ]:scans$Stop[ iii ] ##posn + ii
    chr <- levels( scans$Seq )[ scans$Seq[ iii ] ] ## gene
    m.tmp[[ chr ]][ inds ] <- TRUE ##scans$posns[j]:(scans$posns[j]+wi-1),1] = 1
  }
  cat( i, length( motifs ), motifs[ i ], nrow( scans ), "\n" ) ##sum( m.tmp ), "\n" )
  lapply( m.tmp, which ) ##which( m.tmp ) ##Matrix( m.tmp )
}, mc.preschedule=F )
rm( m.tmp, p.scans ); gc()

if ( is.list( m[[ 1 ]] ) ) {
  for ( i in 1:length( m ) ) {
    x <- m[[ i ]]
    if ( length( unlist( x ) ) <= 0 ) next
    for ( ii in 1:length( x ) ) x[[ ii ]] <- x[[ ii ]] + nc.cumsum[ ii ]
    m[[ i ]] <- unlist( x ); names( m[[ i ]] ) <- NULL
  }
}

#if ( exists( 'DO_STOP' ) && DO_STOP ) stop()
} ## if not exists('m')

#############

e$dlf( 'RegulonDB/RegDB_PSSMSet.txt', 'http://regulondb.ccg.unam.mx/data/PSSMSet.txt' )
tmp <- readLines( 'RegulonDB/RegDB_PSSMSet.txt' )
tfs <- grep( '^Transcription Factor Name', tmp, perl=T )
bs <- as.integer(sapply( strsplit( tmp[ grep( '^Total of binding sites', tmp, perl=T ) ], ' ' ), '[', 5 ) )
names( bs ) <- tfs
pssm.start <- grep( '^a\t\\|\t', tmp, perl=T )
pssms <- list()
for ( i in 1:length(tfs) ) {
  tf <- strsplit( tmp[tfs[i]], ' ' )[[ 1 ]][ 4 ]
  pssm <- do.call(rbind,strsplit(tmp[pssm.start[i]+(0:3)],'\t'))
  cnames <- toupper(pssm[,1])
  pssm <- t( apply( pssm[,-(1:2)], 2, function(i) {i=as.numeric(i);(i+0.01)/sum(i+0.01)} ) )
  colnames( pssm ) <- cnames
  pssms[[tf]] <- pssm
}
names( bs ) <- names( pssms )
cat( "GOT", length(pssms), "ECO PSSMS!!!\n" )
cat( sum(bs>3), "PSSMS BASED ON > 3 SITES\n" ) ## 86
cat( sum(bs>4), "PSSMS BASED ON > 4 SITES\n" ) ## 72

all.pssms.to.mast.file <- function( pssms, seq.type=names(e$mot.weights)[1] ) {
  meme.let <- c( "A", "C", "G", "T" )
  lines <- c( "MEME version 3.0", "", "ALPHABET= ACGT", "", "strands: + -", "",
             "Background letter frequencies (from dataset with add-one prior applied):" )
  lines <- c( lines, paste( names( unlist( e$genome.info$bg.list[[ seq.type ]][ meme.let ] ) ),
                       sprintf( "%.3f", unlist( e$genome.info$bg.list[[ seq.type ]][ meme.let ] ) ), collapse=" " ) )

  cluster.motif.lines <- function( pssm, tf ) { ## Write out a single bicluster's (usually 2) motifs to MEME format
    lines <- character()
    lines <- c( lines, "",
               sprintf( "MOTIF %s", tf ),
               sprintf( "BL   MOTIF %s width=0 seqs=0", tf ),
               sprintf( "letter-probability matrix: alength= 4 w= %d nsites= %d E= %.3e", nrow( pssm ), 10, 1 ) )
    lines <- c( lines, apply( pssm, 1, function( i )
                             sprintf( "%5.3f %5.3f %5.3f %5.3f", i[ 1 ], i[ 2 ], i[ 3 ], i[ 4 ] ) ) )
    lines
  }

  lines.t <- c( lines, do.call( c, lapply( names( pssms ), function( tf ) cluster.motif.lines( pssms[[tf]], tf ) ) ) )
  tfile <- e$my.tempfile( "tomtom_t_", )
  cat( lines.t, file=tfile, sep="\n" ) ## Write out all motifs
  tfile
}

try( load( sprintf( 'fimo_eco_%s.RData', as.character( p.cutoff ) ) ) )
if ( ! exists( 'fimo.eco' ) ) {
  seqs.file <- e$my.tempfile( './fimo_out_eco' )
  writeLines( paste( paste( ">", names( e$genome.info$genome.seqs ), sep="" ), e$genome.info$genome.seqs, sep="\n" ),
             con=seqs.file )
  mots.file <- all.pssms.to.mast.file( pssms ) ## NOTE THIS IS FOR ALL 86 PSSMS - SHOULD FILTER DOWN TO THE 72 GOOD ONES!
  out.file <- e$my.tempfile()

  ## Assume using customized version of meme_4.3.0 with MAX_MOTIFS in motif.h changed to 24000
  fimo.eco <- system( sprintf( './progs/fimo --max-stored-scores 9999999999 --text --verbosity 2 %s %s |bzip2 -c >%s.bz2',
                         mots.file, seqs.file, out.file ), intern=T )
  fimo.eco <- system( sprintf( 'bunzip2 -c %s.bz2 | awk \'($6<=%s){print}\'', out.file, as.character( p.cutoff ) ), intern=T )
  fimo.eco <- as.data.frame( do.call( rbind, strsplit( fimo.eco, '\t' ) ) )
  colnames( fimo.eco ) <- strsplit( system( sprintf( 'bunzip2 -c %s.bz2 | head -1', out.file ), intern=T ), '\t' )[[ 1 ]]
  fimo.eco$Start <- as.integer( as.character( fimo.eco$Start ) )
  fimo.eco$Stop <- as.integer( as.character( fimo.eco$Stop ) )
  fimo.eco$`Log-odds` <- as.numeric( as.character( fimo.eco$`Log-odds` ) )
  fimo.eco$`p-value` <- as.numeric( as.character( fimo.eco$`p-value` ) )
  fimo.eco$Strand <- substr( as.character( fimo.eco$Motif ), 1, 1 )
  fimo.eco$Motif <- substring( fimo.eco$Motif, 2 )
  fimo.eco <- as.data.table( fimo.eco )
  setkey( fimo.eco, Motif, Seq, Start )
  save( fimo.eco, file=sprintf( 'fimo_eco_%s.RData', as.character( p.cutoff ) ) )
}

try( load( sprintf( 'm_eco_%s.RData', as.character( p.cutoff ) ) ) )
if ( ! exists( 'm.eco' ) ) {
  p.scans <- subset( fimo.eco, `p-value` <= p.cutoff )

  m.tmp <- lapply( nc, function( i ) rep( FALSE, length=i ) ); names( m.tmp ) <- names( nc )
  m.eco <- mclapply( names( pssms ), function( tf ) {
    scans <- p.scans[ tf ] ##c( bi, bi ), c( mo, -mo ) ) ]
    if ( nrow(scans) <= 0 ) return( lapply( m.tmp, function( i ) integer(0) ) ) ##Matrix( m.tmp ) )
    scans <- scans[ ! is.na( scans$Start ), ] ##posns ), ]
    if ( nrow(scans) <= 0 ) return( lapply( m.tmp, function( i ) integer(0) ) ) ##Matrix( m.tmp ) )
    for ( iii in 1:nrow( scans ) ) {
      inds <- scans$Start[ iii ]:scans$Stop[ iii ] ##posn + ii
      chr <- levels( scans$Seq )[ scans$Seq[ iii ] ] ## gene
      m.tmp[[ chr ]][ inds ] <- TRUE ##scans$posns[j]:(scans$posns[j]+wi-1),1] = 1
    }
    cat( tf, which( names( pssms ) == tf ), length( pssms ), nrow( scans ), "\n" ) ##sum( m.tmp ), "\n" )
    lapply( m.tmp, which ) ##which( m.tmp ) ##Matrix( m.tmp )
  }, mc.preschedule=F )
  rm( m.tmp, p.scans ); gc()
    
  names( m.eco ) <- names( pssms )
  save( m.eco, file=sprintf( 'm_eco_%s.RData', as.character( p.cutoff ) ) )
}

if ( is.list( m.eco[[ 1 ]] ) ) {
  for ( i in 1:length( m.eco ) ) {
    x <- m.eco[[ i ]]
    if ( length( unlist( x ) ) <= 0 ) next
    for ( ii in 1:length( x ) ) x[[ ii ]] <- x[[ ii ]] + nc.cumsum[ ii ]
    m.eco[[ i ]] <- unlist( x ); names( m.eco[[ i ]] ) <- NULL
  }
}

tmp <- mclapply( names( m.eco ), function( i ) {
  print( i )
  x <- m.eco[[ i ]]
  if ( length( unlist( x ) ) <= 0 ) return()
  tout <- rep( 1, length( m ) )
  for ( j in 1:length( m ) ) {
    y <- m[[ j ]]
    if ( length( y ) <= 0 ) next
    tmp <- x %in% y
    if ( ! any( tmp ) ) next
    tmp <- ( sum( ! tmp ) + sum( ! ( y %in% x ) ) ) / length( unique( c( x, y ) ) ) ##sum( x != y ) / sum( x | y )
    tout[ j ] <- tmp
  }
  tmp <- which( tout < 1 )
  if ( length( tmp ) <= 0 ) return()
  ##write.table( cbind( i, tmp, out[tmp] ), quote=F, sep='\t', row.names=F, col.names=F, file=bzfile( fname ) )
  data.frame( i, motifs[tmp], tout[tmp] )
}, mc.preschedule=F )
eco.hits <- do.call( rbind, tmp )
rm( tmp )
colnames( eco.hits ) <- c( 'eco.tf', 'motif', 'distance' )

## Now try comparing motif.cluster hits to actual locations in RegulonDB
if ( ! exists( 'm.eco2' ) ) {
  tmp <- readLines( 'RegulonDB/BindingSiteSet.txt' )
  if ( FALSE ) {
    tmp <- readLines( 'RegulonDB/BindingSitePredictionSet.txt' )
    skip <- grep( '#   (7) Method', tmp, fixed=T ) + 2
  }
  skip <- grep( 'Evidence that supports the existence', tmp, fixed=T ) + 2
  rm( tmp )
  tab <- read.delim( 'RegulonDB/BindingSiteSet.txt', skip=skip, head=F )
  tab <- subset( tab, V4 != 0 & V5 != 0 )
  n.tot2 <- sort( table( as.character( tab$V2 ) ), decreasing=T )
  n.tot2 <- n.tot2[ n.tot2 >= 3 ]
  tab <- subset( tab, V2 %in% names( n.tot2 ) )
  colnames( tab )[ c( 2, 4, 5 ) ] <- c( 'Motif', 'Start', 'Stop' )
  tab$Strand <- ifelse( tab$V6 == 'forward', '+', '-' )
  tab$Seq <- as.factor( rep( 'NC_000913', nrow( tab ) ) )

  m.tmp <- lapply( nc, function( i ) rep( FALSE, length=i ) ); names( m.tmp ) <- names( nc )
  m.eco2 <- mclapply( names( n.tot2 ), function( tf ) {
    scans <- subset( tab, Motif == tf )
    if ( nrow(scans) <= 0 ) return( lapply( m.tmp, function( i ) integer(0) ) ) ##Matrix( m.tmp ) )
    for ( iii in 1:nrow( scans ) ) {
      inds <- scans$Start[ iii ]:scans$Stop[ iii ] ##posn + ii
      chr <- levels( scans$Seq )[ scans$Seq[ iii ] ] ## gene
      m.tmp[[ chr ]][ inds ] <- TRUE ##scans$posns[j]:(scans$posns[j]+wi-1),1] = 1
    }
    cat( tf, which( names( n.tot2 ) == tf ), length( n.tot2 ), nrow( scans ), "\n" ) ##sum( m.tmp ), "\n" )
    lapply( m.tmp, which ) ##which( m.tmp ) ##Matrix( m.tmp )
  }, mc.preschedule=F )
  names( m.eco2 ) <- names( n.tot2 )
  rm( m.tmp )
}

tmp <- mclapply( names( m.eco2 ), function( i ) {
  print( i )
  x <- m.eco2[[ i ]][[ 2 ]]
  if ( length( x ) <= 0 ) return()
  tout <- rep( 1, length( m ) )
  for ( j in 1:length( m ) ) {
    y <- m[[ j ]] ##[[ 1 ]]
    if ( length( y ) <= 0 ) next
    tmp <- x %in% y
    if ( ! any( tmp ) ) next
    tmp <- ( sum( ! tmp ) + sum( ! ( y %in% x ) ) ) / length( unique( c( x, y ) ) ) ##sum( x != y ) / sum( x | y )
    tout[ j ] <- tmp
  }
  tmp <- which( tout < 1 )
  if ( length( tmp ) <= 0 ) return()
  ##write.table( cbind( i, tmp, out[tmp] ), quote=F, sep='\t', row.names=F, col.names=F, file=bzfile( fname ) )
  data.frame( i, motifs[tmp], tout[tmp] )
}, mc.preschedule=F )
eco.hits2 <- do.call( rbind, tmp )
rm( tmp )
colnames( eco.hits2 ) <- c( 'eco.tf', 'motif', 'distance' )

## Get e.values (need to filter out really high e-val motifs, somehow)
ms <- e$meme.scores[[ seq.type ]]
e.vals <- do.call( rbind, lapply( 1:e$k.clust, function( i ) {
  out <- c( NA, NA, NA )
  if ( is.null( ms[[ i ]] ) || ms[[ i ]] == '' ) return( out )
  mm <- ms[[ i ]]$meme.out
  for ( j in 1:length( mm ) ) if ( ! is.null( mm[[ j ]]$e.value ) ) out[ j ] <- mm[[ j ]]$e.value
  out
} ) )
if ( all( is.na( e.vals[ ,3 ] ) ) ) e.vals <- e.vals[ ,1:2 ]
rm( ms )

## Get combined matches (vs. locations in regdb AND vs. fimo scans of regdb pssms) against 402 best (in terms of e-value) motifs from the single run
max.motifs <- 402 ## Maximum # of motifs per clustering to allow (same as # of motif clusters in eco ensemble)
eval.cutoff <- Inf
if ( sum( ! is.na(e.vals) ) > max.motifs ) eval.cutoff <- sort( e.vals[!is.na(e.vals)] )[ max.motifs ]
tmp <- subset( eco.hits2, distance <= 0.99 ) ## we used 0.99 for the cutoff for the ensemble analysis, but that also had a p-value and frac cutoff so use 0.98 here
tmp <- tmp[ order( tmp$distance ), ]
b.m <- t( apply( sapply( strsplit( as.character( tmp$motif ), '_' ), '[', 2:3 ), 2, as.integer ) )
tmp <- tmp[ e.vals[ b.m ] <= eval.cutoff, ] # use 100 as e-value cutoff (we did not use this for ensemble analysis; what if we did?)
tmp <- subset( tmp, ! duplicated( tmp$motif ) )
print( length( unique( tmp$eco.tf ) ) )


