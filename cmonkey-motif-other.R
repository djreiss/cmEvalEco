
all.dna.seqs <- function( l, lett=c( "G", "A", "T", "C" ), as.matrix=F ) {
  n.lett <- length( lett )
  out <- sapply( 1:l, function( ll ) rep( as.vector( sapply( lett, function( i ) rep( i, n.lett^( ll - 1 ) ) ) ),
                                         n.lett^( l - ll ) ) )
  if ( as.matrix ) return( out )
  apply( out, 1, paste, collapse="" )
}

pssm.motif.lines <- function( pssm, id, e.value=1, header=T, seq.type="upstream weeder" ) {
  meme.let <- c( "A", "C", "G", "T" )
  if ( missing( id ) ) id <- paste( pssm, collapse="" )
  lines <- character()
  if ( header ) {
    lines <- ##c( "MEME version 3.0", "",
               "ALPHABET= ACGT"##, "", "strands: + -", "",
               ##"Background letter frequencies (from dataset with add-one prior applied):" )
               ##)
    lines <- c( lines, paste( names( unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ),
                        sprintf( "%.3f", unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ), collapse=" " ) )
  }
    
  if ( is.null( colnames( pssm ) ) ) colnames( pssm ) <- col.let
  pssm <- pssm[ ,meme.let ] + max( pssm, na.rm=T ) / 100 ## Pseudocount?
  for ( i in 1:nrow( pssm ) ) pssm[ i, ] <- pssm[ i, ] / sum( pssm[ i, ], na.rm=T )
  idd <- gsub( "[_/]", ".", id )
  ##lines <- c( lines, "", sprintf( "MOTIF %s", idd ), sprintf( "BL   MOTIF %s width=0 seqs=0", idd ) )
  ##lines <- c( lines, sprintf( "letter-probability matrix: alength= 4 w= %d nsites= %d E= %.3e", nrow( pssm ),
  ##                           20, e.value ) )
  lines <- c( lines, sprintf( "log-odds matrix: alength= 4 w= %d", nrow( pssm ) ) )
  for ( j in 1:nrow( pssm ) ) lines <- c( lines, paste( sprintf( "%5.3f", log2( pssm[ j, ] ) ),
                                                       collapse=" ", sep=" " ) )
  lines
}

## Write out a single bicluster's (usually 2) motifs to MEME format
cluster.meme.motif.lines <- function( k, seq.type=names( meme.scores )[ 1 ], logodds=F ) { 
  meme.let <- c( "A", "C", "G", "T" )
  lines <- c( "ALPHABET= ACGT", "", "strands: + -", "",
             "Background letter frequencies (from dataset with add-one prior applied):" )
  lines <- c( lines, paste( names( unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ),
                       sprintf( "%.3f", unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ), collapse=" " ) )
  memeOut <- meme.scores[[ seq.type ]][[ k ]]
  if ( is.null( memeOut ) || ! is.list( memeOut ) || is.null( memeOut$meme.out ) ) return( lines )
  memeOut <- memeOut$meme.out
  for ( i in 1:length( memeOut ) ) {
    pssm <- memeOut[[ i ]]$pssm ##; colnames( pssm ) <- meme.let ##col.let; pssm <- pssm[ ,meme.let ]
    mat.type <- "letter-probability matrix"
    if ( logodds ) mat.type <- "log-odds matrix"
    lines <- c( lines, "",
               sprintf( "MOTIF bic_%03d_%02d", k, i ),
               sprintf( "BL   MOTIF bic_%03d_%02d width=0 seqs=0", k, i ),
               sprintf( "%s: alength= 4 w= %d nsites= %d E= %.3e", mat.type, nrow( pssm ),
                       memeOut[[ i ]]$sites, memeOut[[ i ]]$e.value ) )
    if ( ! logodds ) lines <- c( lines, apply( pssm, 1, function( i )
                                              sprintf( "%5.3f %5.3f %5.3f %5.3f", i[ 1 ], i[ 2 ], i[ 3 ], i[ 4 ] ) ) )
    else lines <- c( lines, apply( round( log( pssm + 0.01 ) ), 1, function( i )
                                  sprintf( "%6d %6d %6d %6d", i[ 1 ], i[ 2 ], i[ 3 ], i[ 4 ] ) ) )
  }
  lines
}

all.motifs.to.mast.file <- function( ks=1:k.clust, seq.type=names(mot.weights)[1],
                                    e.value.cutoff=100, resid.cutoff=0.8 ) {
  meme.let <- c( "A", "C", "G", "T" )
  lines <- c( "MEME version 3.0", "", "ALPHABET= ACGT", "", "strands: + -", "",
             "Background letter frequencies (from dataset with add-one prior applied):" )
  lines <- c( lines, paste( names( unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ),
                       sprintf( "%.3f", unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ), collapse=" " ) )

  cluster.motif.lines <- function( k ) { ## Write out a single bicluster's (usually 2) motifs to MEME format
    if ( k %% 100 == 0 ) print(k)
    lines <- character()
    memeOut <- meme.scores[[ seq.type ]][[ k ]]
    if ( is.null( memeOut ) || memeOut == "" ) return( lines )
    memeOut <- meme.scores[[ seq.type ]][[ k ]]$meme.out
    if ( is.null( memeOut ) ) return( lines )
    if ( clusterStack[[ k ]]$resid > resid.cutoff ) return( lines )
    ##max.motifs <- max( max.motifs, length( memeOut ) )
    for ( i in 1:length( memeOut ) ) {
      if ( memeOut[[ i ]]$e.value > e.value.cutoff ) next
      pssm <- memeOut[[ i ]]$pssm ##; colnames( pssm ) <- meme.let ##col.let; pssm <- pssm[ ,meme.let ]
      lines <- c( lines, "",
                 sprintf( "MOTIF bic_%03d_%02d_%.3f_%.3e", k, i, clusterStack[[ k ]]$resid, memeOut[[ i ]]$e.value ),
                 sprintf( "BL   MOTIF bic_%03d_%02d_%.3f_%.3e width=0 seqs=0", k, i, clusterStack[[ k ]]$resid, memeOut[[ i ]]$e.value ),
                 sprintf( "letter-probability matrix: alength= 4 w= %d nsites= %d E= %.3e", nrow( pssm ),
                         memeOut[[ i ]]$sites, memeOut[[ i ]]$e.value ) )
      lines <- c( lines, apply( pssm, 1, function( i )
                               sprintf( "%5.3f %5.3f %5.3f %5.3f", i[ 1 ], i[ 2 ], i[ 3 ], i[ 4 ] ) ) )
    }
    lines
  }

  lines.t <- c( lines, do.call( c, lapply( ks, cluster.motif.lines ) ) )
  tfile <- my.tempfile( "tomtom_t_", )
  cat( lines.t, file=tfile, sep="\n" ) ## Write out all motifs
  tfile
}  

## All vs. all comparison of motifs in every cluster using MEME package's tomtom program
## Note that tomtom can only handle up to 5000 motifs at a time!!!
## But this can be changed in the tomtom.c code -- I have done it on pinnacle to go up to 15000 motifs
motif.similarities.tomtom <- function( query=1:k.clust, target=1:k.clust, query.mot=NA, target.mot=NA,
                                      seq.type=names(mot.weights)[1],
                                      e.value.cutoff=100, resid.cutoff=0.8, dist.meth="ed", q.thresh=0.5,
                                      min.overlap=4, q.pseudo=0, t.pseudo=0, min.gene.overlap=NA,
                                      desymmetrize=T, unlink=T, files.only=F, verbose=T, save.touts=F, ... ) {
  meme.let <- c( "A", "C", "G", "T" )
  lines <- c( "MEME version 3.0", "", "ALPHABET= ACGT", "", "strands: + -", "",
             "Background letter frequencies (from dataset with add-one prior applied):" )
  lines <- c( lines, paste( names( unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ),
                       sprintf( "%.3f", unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ), collapse=" " ) )

  query <- query[ ! is.na( query ) ]
  target <- target[ ! is.na( target ) ]
  if ( is.na( query.mot ) ) query.mot <- rep( NA, length( query ) )
  if ( is.na( target.mot ) ) target.mot <- rep( NA, length( target ) )
  ##max.motifs <- 0
  
  cluster.motif.lines <- function( k, mot ) { ## Write out a single bicluster's (usually 2) motifs to MEME format
    lines <- character()
    memeOut <- meme.scores[[ seq.type ]][[ k ]]
    if ( is.null( memeOut ) || memeOut == "" ) return( lines )
    memeOut <- memeOut$meme.out ##meme.scores[[ seq.type ]][[ k ]]$meme.out
    if ( is.null( memeOut ) ) return( lines )
    clust <- clusterStack[[ k ]]
    if ( is.na( clust$resid ) || clust$resid > resid.cutoff ) return( lines )
    ##max.motifs <- max( max.motifs, length( memeOut ) )
    for ( i in 1:length( memeOut ) ) {
      if ( ! is.na( mot ) && length( mot ) > 0 && ! i %in% mot ) next
      mo <- memeOut[[ i ]]
      if ( is.na( mo$e.value ) || mo$e.value > e.value.cutoff ) next
      pssm <- mo$pssm ##; colnames( pssm ) <- meme.let ##col.let; pssm <- pssm[ ,meme.let ]
      lines <- c( lines, "",
                 sprintf( "MOTIF bic_%03d_%02d_%.3f_%.3e", k, i, clust$resid, mo$e.value ),
                 sprintf( "BL   MOTIF bic_%03d_%02d_%.3f_%.3e width=0 seqs=0", k, i, clust$resid, mo$e.value ),
                 sprintf( "letter-probability matrix: alength= 4 w= %d nsites= %d E= %.3e", nrow( pssm ),
                         mo$sites, mo$e.value ) )
      lines <- c( lines, apply( pssm, 1, function( i )
                               sprintf( "%5.3f %5.3f %5.3f %5.3f", i[ 1 ], i[ 2 ], i[ 3 ], i[ 4 ] ) ) )
    }
    lines
  }

  cmd <- "%s/tomtom -verbosity 1 -q-thresh %.3f -dist %s -min-overlap %d -text -query-pseudo %.3f -target-pseudo %.3f -target %s"
  mc <- get.parallel( length( target ) ) ## Parallelize running query motifs against all target motifs
  apply.fun <- mc$apply
  if ( multicore:::isChild() || length( query ) == 1 ) apply.func <- lapply
  ##apply.fun <- lapply
  
  if ( is.na( min.gene.overlap ) ) { ## If no cluster-specific targets, just create target file once for all clusters
    ##if ( files.only == TRUE || ( ! is.logical( files.only ) && ! is.na( files.only ) && ! is.null( files.only ) ) )
    ##  apply.fun <- mclapply
    lines.t <- c( lines, unlist( apply.fun( target, function( k )
                                          cluster.motif.lines( k, target.mot[ target == k ] ) ) ) )
    if ( verbose ) cat( "TARGET MOTIFS:", range( target ), length( grep( "BL   MOTIF", lines.t, fixed=T ) ), "\n" )
    tfile <- my.tempfile( "tomtom_t_", ); cat( lines.t, file=tfile, sep="\n" ) ## Write out all target motifs
    cmd <- sprintf( cmd, progs.dir, q.thresh, dist.meth, min.overlap, q.pseudo, t.pseudo, tfile )
  }

  if ( ! is.na( min.gene.overlap ) ) c.rows <- lapply( 1:k.clust, get.rows )
  if ( files.only == TRUE || ( ! is.logical( files.only ) && ! is.na( files.only ) && ! is.null( files.only ) ) ) ##{
    options( cores=1 )
  if ( save.touts != FALSE ) dir.create( save.touts, recursive=T ) ##"filehash/touts", recursive=T )
  tout <- apply.fun( unique( query ), function( k ) {
    if ( save.touts != FALSE ) { ##min( ks, na.rm=T ) ) ) ) {
      tout.file <- sprintf( "%s/%08d_tout.tsv.bz2", save.touts, k )
      if ( file.exists( tout.file ) ) {
        print( tout.file ) ##sprintf( "filehash/touts/%08d_tout.tsv.bz2", k ) ) ##min( ks, na.rm=T ) ) )
        tout <- try( read.delim( bzfile( tout.file ), sep='\t', head=T, check.names=F ) )
        if ( class( tout ) != 'try-error' ) {
          if ( ! is.infinite( e.value.cutoff ) && ! is.na( e.value.cutoff ) ) {
            e.val1 <- as.numeric( sapply( strsplit( as.character( tout[ ,1 ] ), "_" ), "[", 5 ) )
            e.val2 <- as.numeric( sapply( strsplit( as.character( tout[ ,2 ] ), "_" ), "[", 5 ) )
            tout <- subset( tout, e.val1 <= e.value.cutoff & e.val2 <= e.value.cutoff )
          }
          print( dim( tout ) )
          return( tout )
        }
      }
    }
    lines.q <- c( lines, cluster.motif.lines( k, query.mot[ which( query == k ) ] ) )## ) )
    n.query.motifs <- length( grep( "BL   MOTIF", lines.q, fixed=T ) )
    if ( verbose ) cat( "QUERY MOTIFS:", k, n.query.motifs, length( which( query == k ) ), "\n") ##range( ks ), "\n" )
    if ( n.query.motifs <= 0 ) return( NULL )
    if ( ! is.na( min.gene.overlap ) ) {
      rows <- c.rows[[ k ]] ##unlist( c.rows[ k ] ) ##s ] )
      t.ok <- sapply( target, function( t ) sum( c.rows[[ t ]] %in% rows ) >= min.gene.overlap )
      t.ok[ which( target == k ) ] <- FALSE ##%in% ks ) ] <- FALSE
      lines.t <- c( lines, unlist( lapply( target[ t.ok ], function( kk )
                                             cluster.motif.lines( kk, target.mot[ which( target == kk ) ] ) ) ) )
      n.target.motifs <- length( grep( "BL   MOTIF", lines.t, fixed=T ) )
      if ( n.target.motifs <= 0 ) {
        tout <- NULL
        if ( ! files.only && save.touts != FALSE && ! file.exists( sprintf( "%s/%08d_tout.tsv.bz2", save.touts, k ) ) ) {
          writeLines( '', con=gsub( '.bz2', '', tout.file, fixed=T ) ) ##min( ks, na.rm=T ) ) )
          system( sprintf( "bzip2 -fv %s/%08d_tout.tsv", save.touts, k ) ) ##min( ks, na.rm=T ) ) )
        }
        return( NULL )
      }
      if ( verbose ) cat( "TARGET MOTIFS:", k, n.target.motifs, sum( t.ok ), range( target[ t.ok ] ), "\n" )
      tfile <- my.tempfile( "tomtom_t_", ); cat( lines.t, file=tfile, sep="\n" ) ## Write out all target motifs
      cmd <- sprintf( cmd, progs.dir, q.thresh, dist.meth, min.overlap, q.pseudo, t.pseudo, tfile )
    }
    ##qfile <- my.tempfile( "tomtom_q_" )
    qfile <- paste( gsub( "tomtom_t_", "tomtom_q_", tfile ), "_", k, '_', sep="" )
    cat( lines.q, file=qfile, sep="\n" )
    cmd <- paste( cmd, "-query", qfile )
    if ( files.only == TRUE || ( ! is.logical( files.only ) && ! is.na( files.only ) && ! is.null( files.only ) ) ) {
      out.cmd <- paste( cmd, " > ", qfile, ".out", sep="" )
      if ( is.character( files.only ) ) cat( out.cmd, "\n", file=files.only, append=T )
      if ( verbose ) print( out.cmd )
      return( out.cmd )
    }
    if ( verbose ) print( cmd )
    tout <- system( cmd, intern=T )
    if ( unlink ) {
      if ( ! is.na( min.gene.overlap ) ) unlink( tfile )
      unlink( qfile )
    }
    if ( ! is.null( tout ) ) {
      tout <- do.call( rbind, strsplit( tout, "\t" ) )
      colnames( tout ) <- tout[ 1, ,drop=F ]; tout <- tout[ -1, ,drop=F ]
      tout <- as.data.frame( tout[ tout[ ,1 ] != tout[ ,2 ], ,drop=F ] ) ## < ##de-symmetrize the output )upper-tri)
      cat( "TOUT:", dim( tout ), "\n" )
      if ( nrow( tout ) > 0 && save.touts != FALSE && ! file.exists( sprintf( "%s/%08d_tout.tsv.bz2", save.touts, k ) ) ) {
        write.table( tout, file=sprintf( "%s/%08d_tout.tsv", save.touts, k ), ##min( ks, na.rm=T ) ),
                    quote=F, sep='\t', append=F, row.names=F )
        system( sprintf( "bzip2 -fv %s/%08d_tout.tsv", save.touts, k ) ) ##min( ks, na.rm=T ) ) )
      }
    } else {
      print( "NULL TOUT!" )
      warning( "NULL TOUT!" )
    }
    tout
  } ) ##, mc.preschedule=F )

  if ( ! is.data.frame( tout ) ) cat( "GOT", length( tout ), "touts.\n" )
  if ( files.only == TRUE || ( ! is.logical( files.only ) && ! is.na( files.only ) && ! is.null( files.only ) ) )
    return( tout )
  if ( ! is.data.frame( tout ) ) {
    cat( "NROW:", sum( unlist( sapply( tout, nrow ) ) ), "\n" )
    if ( sum( unlist( sapply( tout, nrow ) ) ) > 0 ) {
      tout <- tout[ which( ! sapply( tout, function(i) is.null(i) || class(i)=='try-error' ) ) ]
      tout <- tout[ which( sapply( tout, nrow ) > 0 ) ]
      ##tout2 <- do.call( rbind, tout )
      cn <- colnames( tout[[ 1 ]] )
      for ( i in 1:length( tout ) ) colnames( tout[[ i ]] ) <- rownames( tout[[ i ]] ) <- NULL
      for ( i in 1:length( tout ) ) for ( j in 1:ncol( tout[[ i ]] ) )
        tout[[ i ]][ ,j ] <- as.character( tout[[ i ]][ ,j ] )
      for ( i in 1:length( tout ) ) tout[[ i ]] <- as.matrix( tout[[ i ]] )
      gc()
      tout <- do.call( rbind, tout )
      colnames( tout ) <- cn
      tout <- as.data.frame( tout )
    } else {
      tout <- tout[[ 1 ]]
    }
  } 

  
  cat( "GOT", nrow( tout ), "motif alignments.\n" )
  if ( is.null( tout ) ) return( NULL )
  colnames( tout ) <- gsub( '^X\\.', '', colnames( tout ) )
  ## Summarize into data frame
  if ( nrow( tout ) > 0 ) {
    q.id <- t( sapply( strsplit( as.character( tout[ ,1 ] ), "_", fixed=T ), '[', 2:5 ) ) ##function( i ) i[ 2:5 ] )
    ##q.res.ev <- do.call( rbind, strsplit( as.character( tout[ ,1 ] ), "_" ), function( i ) as.numeric( i[ 4:5 ] ) )
    t.id <- t( sapply( strsplit( as.character( tout[ ,2 ], fixed=T ), "_" ), '[', 2:5 ) )
    ##t.res.ev <- do.call( rbind, strsplit( as.character( tout[ ,2 ] ), "_" ), function( i ) as.numeric( i[ 4:5 ] ) )
    tout2 <- data.frame( biclust1=as.integer( q.id[ ,1 ] ), motif1=as.integer( q.id[ ,2 ] ),
                        resid1=as.numeric( q.id[ ,3 ] ), e.value1=as.numeric( q.id[ ,4 ] ),
                        biclust2=as.integer( t.id[ ,1 ] ), motif2=as.integer( t.id[ ,2 ] ),
                        resid2=as.numeric( t.id[ ,3 ] ), e.value2=as.numeric( t.id[ ,4 ] ),
                        offset=as.integer( as.character( tout[[ 'Optimal offset' ]] ) ),
                        p.value=as.numeric( as.character( tout[[ 'p-value' ]] ) ),
                        q.value=as.numeric( as.character( tout[[ 'q-value' ]] ) ),
                        overlap=as.integer( as.character( tout[[ 'Overlap' ]] ) ),
                        consensus1=as.factor( as.character( tout[[ 'Query consensus' ]] ) ),
                        consensus2=as.factor( ifelse( tout[[ 'Orientation' ]] == "-",
                          rev.comp( as.character( tout[[ 'Target consensus' ]] ) ),
                          as.character( tout[[ 'Target consensus' ]] ) ) ),
                        orientation=tout[[ 'Orientation' ]] ) ##,
    tout2 <- subset( tout2, ! ( biclust1 == biclust2 & motif1 == motif2 ) )
    tout2 <- tout2[ order( tout2$q.value, tout2$p.value ), ]
    if ( desymmetrize ) tout2 <- desymmetrize.tomtom.results( tout2 )
  } else {
    tout2 <- data.frame( biclust1=integer(), motif1=integer(), resid1=numeric(), e.value1=numeric(),
                        biclust2=integer(), motif2=integer(), resid2=numeric(), e.value2=numeric(),
                        offset=integer(), p.value=numeric(), q.value=numeric(), overlap=integer(),
                        consensus1=factor(), consensus2=factor(), orientation=factor() )
  }
  rm( tout )  
  ##else tout <- tout2
  if ( exists( cmd ) ) attr( tout2, "tomtom.cmd" ) <- cmd
  tout2
}

## lower tri != upper tri because of asymmetries in background - so lets use the min of a vs b and b vs a:
## do it matrix-y by min-ing the upper- and lower- tri of a p-value matrix
desymmetrize.tomtom.results <- function( tt.out ) {
  mot.names <- unique( c( paste( tt.out$biclust1, tt.out$motif1, sep="_" ),
                         paste( tt.out$biclust2, tt.out$motif2, sep="_" ) ) )
  mot.names.2 <- cbind( paste( tt.out$biclust1, tt.out$motif1, sep="_" ),
                       paste( tt.out$biclust2, tt.out$motif2, sep="_" ) )
  mot.names.2a <- cbind( mot.names.2[ ,2 ], mot.names.2[ ,1 ] )
  tmp <- matrix( NA, nrow=length( mot.names ), ncol=length( mot.names ) )
  rownames( tmp ) <- colnames( tmp ) <- mot.names
  lt <- lower.tri( tmp )
  tmp[ mot.names.2 ] <- tt.out$p.value
  tmp2 <- cbind( tmp[ lt ], t( tmp )[ lt ] )
  tmp2[ is.na( tmp2 ) ] <- Inf
  tmp[ lt ] <- apply( tmp2, 1, min, na.rm=T )
  tmp[ is.infinite( tmp ) ] <- NA
  rm( tmp2 )
  tmp.good <- ! is.na( tmp[ lt ] )

  tmp2 <- tmp * NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$biclust1
  out <- data.frame( biclust1=tmp2[ lt ][ tmp.good ] )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$motif1
  out <- cbind( out, data.frame( motif1=tmp2[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$resid1
  out <- cbind( out, data.frame( resid1=tmp2[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$e.value1
  out <- cbind( out, data.frame( e.value1=tmp2[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$biclust2
  out <- cbind( out, data.frame( biclust2=tmp2[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$motif2
  out <- cbind( out, data.frame( motif2=tmp2[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$resid2
  out <- cbind( out, data.frame( resid2=tmp2[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$e.value2
  out <- cbind( out, data.frame( e.value2=tmp2[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$offset
  out <- cbind( out, data.frame( offset=tmp2[ lt ][ tmp.good ] ) )
  out <- cbind( out, data.frame( p.value=tmp[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$q.value
  out <- cbind( out, data.frame( q.value=tmp2[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$overlap
  out <- cbind( out, data.frame( overlap=tmp2[ lt ][ tmp.good ] ) )
  ##if ( is.factor( tt.out$orientation ) ) {
  ##  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- as.numeric( tt.out$orientation )
  ##  out <- cbind( out, data.frame( orientation=tmp2[ lt ][ tmp.good ] ) )
  ##}
  
  rm( tmp2 )
  if ( "consensus1" %in% colnames( tt.out ) || "consensus2" %in% colnames( tt.out ) ||
      is.character( tt.out$orientation ) ) {
    tmp.c <- matrix( "", nrow=length( mot.names ), ncol=length( mot.names ) )
    rownames( tmp.c ) <- colnames( tmp.c ) <- mot.names
    if ( "consensus1" %in% colnames( tt.out ) ) {
      tmp.c[ mot.names.2 ] <- tmp.c[ mot.names.2a ] <- as.character( tt.out$consensus1 )
      out <- cbind( out, data.frame( consensus1=tmp.c[ lt ][ tmp.good ] ) )
    }
    if ( "consensus2" %in% colnames( tt.out ) ) {
      tmp.c[,] <- ""; tmp.c[ mot.names.2 ] <- tmp.c[ mot.names.2a ] <- as.character( tt.out$consensus2 )
      out <- cbind( out, data.frame( consensus2=tmp.c[ lt ][ tmp.good ] ) )
    }
    if ( is.factor( tt.out$orientation ) ) tt.out$orientation <- as.character( tt.out$orientation )
    if ( is.character( tt.out$orientation ) ) {
      tmp.c[,] <- ""; tmp.c[ mot.names.2 ] <- tmp.c[ mot.names.2a ] <- as.character( tt.out$orientation )
      out <- cbind( out, data.frame( orientation=tmp.c[ lt ][ tmp.good ] ) )
    }
    rm( tmp.c )
  }
  
  out <- subset( out, ! is.na( p.value ) )
  out <- out[ order( out$p.value, out$q.value ), ]
  out
}
