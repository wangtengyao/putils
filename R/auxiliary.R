############### Auxiliary functions ###############

#' Print line
#' @param ... R objects printable using cat()
#' @details This is a wrapper function that prints R objects using cat() with
#' no space separation and end with a newline character
#' @export
println <- function(...){
  .Internal(cat(c(list(...), '\n'), file=stdout(), sep='', fill=FALSE, labels=NULL, append=FALSE))
}

#' Write Latex Table
#' @param d if specified, all floats will be rounded to d decimal places
#' @param s if specified (and d is unspecified), then all floats will be rounded
#' to s signif
#' @param no.rounding column indices not to apply rounding
#' @param file filename for output, default to screen output
#' @examples
#' df <- data.frame(name=c('Alpha', 'Beta', 'Gamma', 'Delta'),
#'                  size=c(100L,200L,300L,400L), score=c(23.091,19.978,1119.229, 0.03089))
#' write.latextable(df, s=3)
#' @export
write.latextable <- function(x, d=NA, s=NA, no.rounding=numeric(), file=''){
  float_rounding <- function(a){
    if (!is.na(d)) return(dp(a, digits=d))
    if (!is.na(s)) return(sf(a, digits=s))
    return(a)
  }
  df <- x
  # preprocess columns of df by type
  for (i in seq_along(df)){
    if (class(df[[i]])=='integer'){
      df[[i]] <- paste0('$', as.character(df[[i]]), '$')
    } else if (class(df[[i]])=='numeric' && !(i %in% no.rounding)){
      df[[i]] <- paste0('$', sapply(df[[i]], float_rounding), '$')
    } else {
      df[[i]] <- as.character(df[[i]])
    }
  }
  cat('\\begin{tabular}{', rep('c', ncol(df)), '}\n\\hline\\hline\n', sep='', file=file)
  cat(paste0(colnames(df), collapse=' & '), '\\\\\n\\hline\n', sep='', file=file, append=TRUE)
  write.table(df, file=file, append=TRUE, quote=FALSE, sep=' & ', eol='\\\\\n',
              na=' ', row.names=FALSE, col.names=FALSE)
  cat('\\hline\\hline\n\\end{tabular}\n', file=file, append=TRUE)
}


#' Print percentage
#' @param ind a vector of for loop interator
#' @param tot a vector of for loop lengths
#' @return on screen output of percentage
#' @export
printPercentage <- function (ind, tot){
    ind <- as.vector(ind); tot <- as.vector(tot)
    if ((length(tot) > 1) & (length(ind) == 1)) {ind <- match(ind, tot); tot <- length(tot)}
    len <- length(ind)
    contrib <- rep(1,len)
    if (len > 1) {
        for (i in (len-1):1) contrib[i] <- contrib[i+1] * tot[i+1]
    }
    grand_tot <- contrib[1] * tot[1]
    count <- (sum(contrib * (ind - 1)) + 1)
    out <- ""
    if (sum(ind-1)>0) out <- paste0(rep("\b", nchar(round((count-1)/grand_tot * 100))+1), collapse = "")
    out <- paste0(out, round(count/grand_tot*100), "%")
    if (identical(ind, tot)) out <- paste0(out, '\n')
    cat(out)
    return(NULL)
}

#' Visualise a matrix X
#' @param X a matrix
#' @param aspect.ratio if automatic, it will be calculated automatically to fit screen, otherwise, the actual dimension of the matrix will be used.
#' @param axes whether to display axes
#' @param frame.plot whether to draw a frame
#' @return a color plot of matrix value magnitude
#' @export
visualise <- function(X, aspect.ratio = c('automatic', 'actual'), axes = FALSE, frame.plot = FALSE){
    aspect.ratio = match.arg(aspect.ratio)
    n = dim(X)[1]; p = dim(X)[2]
    if (aspect.ratio == 'actual') {
        image(t(X[n:1,]),asp=n/p, axes = axes, frame.plot = frame.plot)
    }
    else {
        image(t(X[n:1,]), axes = axes, frame.plot = frame.plot)
    }
}

#' Show snippet of a large vector/matrix
#' @param A a vector, matrix or array
#' @param nrow number of rows to show
#' @param ncol number of columns to show, ignored for vectors
#' @details Show the first nrow entries of a vector, the first nrow x ncol
#' submatrix of a matrix. If A is an array, then randomly sample the third to
#' the last indices and show the first nrow x ncol entries in that frame.
#' @export
snippet <- function(A, nrow=5, ncol=nrow){
  if (is.vector(A)){
    cat('Vector of length ', length(A), ', with leading entries:\n', sep='')
    print(A[seq_len(min(length(A), nrow))])
  } else if (is.matrix(A)) {
    cat('Matrix with shape (', paste(as.character(dim(A)), collapse=', '),
        '), with leading entries:\n', sep='')
    print(A[seq_len(min(nrow, nrow(A))), seq_len(min(ncol, ncol(A)))])
  } else if (is.array(A)) {
    dims <- dim(A); d <- length(dims);
    shape <- paste(as.character(dim(A)), collapse=', ')
    if (d == 1){
      cat('1-d array of length ', dims, ', with leading entries:\n', sep='')
      print(A[seq_len(min(length(A), nrow))])
    } else if (d == 2){
      cat('2-d array with shape (', shape, '), with leading entries:\n', sep='')
      print(A[seq_len(min(nrow, nrow(A))), seq_len(min(ncol, ncol(A)))])
    } else {
      frames <- rep(0, d-2); starting_index <- 0
      for (i in seq_len(d-2)){
        frames[d-1-i] <- sample(dims[d+1-i], 1)
        starting_index <- starting_index + prod(head(dims, d-i)) * (frames[d-1-i] - 1)
      }
      cat(d, '-d array with shape (', shape, '), with leading entries in frame [:, :, ',
          paste(as.character(frames), collapse=', '), ']:\n', sep='')
      M <- matrix(A[starting_index + seq_len(dims[1]*dims[2])], dims[1], dims[2])
      print(M[seq_len(min(nrow, nrow(M))), seq_len(min(ncol, ncol(M)))])
    }
  } else {
    stop('A need to be a vector or a matrix or an array.')
  }
}

#' Find the location of first TRUE value in a boolean vector
#' @param v a logical vector
#' @return an integer denotating the location, return NA if not found.
#' @export
find.first <- function(v){
  match(TRUE, v, nomatch = NA)
}

#' Find the location of final TRUE value in a boolean vector
#' @param v a logical vector
#' @return an integer denotating the location, return NA if not found.
#' @export
find.last <- function(v){
  n <- length(v)
  n + 1L - match(TRUE, rev(v), nomatch = NA)
}

#' display signif of exponentiated number nicely
#' @details significant figure computed after subtracting 1. keep trailing zeros, not use scientific notation
#' @param x a real number
#' @param digits positive integer, number of significant figures
#' @return a string
#' @export
sf_exp <- function(x, digits){
  as.character(as.numeric(sf(x-1,digits))+1)
}

#' display signif nicely
#' @details keep trailing zeros, not use scientific notation
#' @param x a real number
#' @param digits number of significant figures to keep
#' @return a string
#' @export
sf <- function(x, digits){
  str <- formatC(signif(x, digits=digits), digits=digits, format="fg", flag="#")
  if (suffix(str, 1)=='.') str <- prefix(str, nchar(str) - 1)
  return(str)
}

#' display decimal places nicely
#' @details keep trailing zeros
#' @param x a real number
#' @param digits number of decimal places to keep
#' @return a string
#' @export
dp <- function(x, digits){
  digits <- floor(log10(x)) + 1 + digits
  str <- formatC(round(x, digits), digits=digits, format="fg", flag = "#")
  if (suffix(str, 1)=='.') str <- prefix(str, nchar(str) - 1)
  return(str)
}

#' Simulation parameter data frame generation
#' @description  create a dataframe of all possible parameter combinations in lexicographic order (if tags are supplied, use tag for column names)
#' @param ... each argument should be of the form of tag = vector, meaning the variable named 'tag' takes values in 'vector'.
#' @details A sample usage is sim.params(tag1 = vec1, tag2 = vec2, tag3 = vec3).
#' @export
sim.params <- function(...){
  x <- list(...)
  n <- length(x)
  vnames <- names(x); no.vn <- !nzchar(vnames)
  vnames[no.vn] <- paste0('Var', seq_len(n))[no.vn]
  df <- expand.grid(rev(x))[,rev(seq_len(n))]
  colnames(df) <- vnames
  return(df)
}

#' Show parameter values
#' @description Print out parameters in a vector in a nice format
#' @param x a vector of parameters
#' @export
show.params <- function(...) {
  names <- as.list(substitute(list(...)))[-1L]
  vals <- list(...)
  paste(paste0(names, ' = ', vals), collapse=', ')
}


#' Multiple assignment
#' @description assign multiple items in a list on RHS to multiple items in a list on LHS
#' @details A sample usage is  \code{bunch(a,b,c) %=% list('hello', 123, list('apple', 'orange'))}, or \code{bunch(a,b,c) %=% 1:3}
#' @param l left side list, enclosed by the \code{bunch} function
#' @param r right side list
#' @export
'%=%' <- function(l, r) UseMethod('%=%')  # Generic form

#' Binary Operator
#' @description method for lbunch
#' @export
'%=%.lbunch' = function(l, r) {
  Envir = as.environment(-1)
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  for (i in seq_along(l)) {
    do.call('<-', list(l[[i]], r[[i]]), envir=Envir)
  }
}

# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)

  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin

  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}

#' Grouping the left hand side in multiple assignment
#' @description bunch multiple items together for multiple assignment
#' @param ... variables to be bunched
#' @return a list of variable names
#' @export
bunch = function(...) {
  List <- as.list(substitute(list(...)))[-1L]
  class(List) <- 'lbunch'
  return(List)
}


#' Change all NA values in v to a
#' @param v a vector
#' @param a target value
#' @return updated vector
#' @export
setNA <- function(v, a){
  v[is.na(v)] <- a;
  return(v)
}

#' matplotlib_palette
#' return the first n palette colours
#' @param n number of colours
#' @param scheme palette scheme, one of 'default', 'bright' and 'rainbow'
#' @param visualise if TRUE, a barplot of all colours will be shown
#' @return a vector of hexadecimal colours
#' @export
matplotlib_palette <- function(n=0, scheme='default', visualise=FALSE){
  default_palette <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
                       "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
  bright_palette <- c("#0e4897","#17813f","#1d99b4","#1f9ee8","#25ca7a","#471c7c",
                      "#68c7ed","#6d4e98","#73af38","#7f1273","#9e1653","#ab0077",
                      "#b01426","#b1b2b4","#c1d430","#cc0b24","#e10064","#e12653",
                      "#e34e9d","#e46b07","#fbee29","#fcc125")
  rainbow_palette <- c("#BF4D4D","#BF864D","#BFBF4D","#86BF4D","#4DBF4D",
                       "#4DBF86","#4DBFBF","#4D86BF","#4D4DBF","#864DBF",
                       "#BF4DBF","#BF4D86")
  full_palette <- switch(scheme,
                         'default' = default_palette,
                         'bright' = bright_palette,
                         'rainbow' = rainbow_palette)

  if (n == 0) {
    ret <- full_palette
  } else {
    reps <- ceiling(n / length(full_palette))
    ret = character()
    for (rep in 1:reps){
      mod_color <- unname(sapply(full_palette, function(c)mix_color(c, '#ffffff', (rep-1)/reps)))
      if (rep==reps) mod_color <- head(mod_color, n - length(full_palette)*(reps-1))
      ret <- c(ret, mod_color)
    }
  }

  if (visualise) barplot(rep(1, length(ret)), col=ret, axes=F, border=F,
                         names.arg=seq_along(ret), cex.names=0.8)
  return(ret)
}

#' Color mixing
#' @param color1 color 1 in hexadecimal
#' @param color2 color 2 in hexadecimal
#' @param lambda interpolation, if equal to 0, return color1, if 1 return color2
#' @export
mix_color <- function(color1, color2, lambda=0.5){
  R1 <- strtoi(paste0('0x', substr(color1, 2, 3)))
  G1 <- strtoi(paste0('0x', substr(color1, 4, 5)))
  B1 <- strtoi(paste0('0x', substr(color1, 6, 7)))
  R2 <- strtoi(paste0('0x', substr(color2, 2, 3)))
  G2 <- strtoi(paste0('0x', substr(color2, 4, 5)))
  B2 <- strtoi(paste0('0x', substr(color2, 6, 7)))
  R <- round(R1 * (1-lambda) + R2 * lambda, 0)
  G <- round(G1 * (1-lambda) + G2 * lambda, 0)
  B <- round(B1 * (1-lambda) + B2 * lambda, 0)
  return(rgb(R, G, B, maxColorValue = 255))
}


#' Set minus: remove elements of small set from large set
#' @param large the large set
#' @param small the small set
#' @return large set with elements of small set removed
`%setminus%` <- function(large, small){
  loc <- match(small, large)
  loc <- loc[!is.na(loc)]
  return(large[-loc])
}


#' Binary operators, add/subtract/multiply/divide a vector to a matrix row by
#' row. i.e. each row of the matrix is added/subtracted/etc by the same vector
#' @param x matrix or vector
#' @param y matrix or vector (one of x and y needs to be a vector of length
#' equal to the number of columns of the other)
#' @name sweep_arithmetic
NULL

#' @rdname sweep_arithmetic
#' @export
`%_+_%` <- function(x, y) {
  if (is.matrix(x) && !is.matrix(y)) return(t(t(x) + y))
  if (is.matrix(y) && !is.matrix(x)) return(t(x + t(y)))
  stop('two arguments must contain exactly one matrix and one vector')
}
#' @rdname sweep_arithmetic
#' @export
`%_-_%` <- function(x, y) {
  if (is.matrix(x) && !is.matrix(y)) return(t(t(x) - y))
  if (is.matrix(y) && !is.matrix(x)) return(t(x - t(y)))
  stop('two arguments must contain exactly one matrix and one vector')
}
#' @rdname sweep_arithmetic
#' @export
`%_*_%` <- function(x, y) {
  if (is.matrix(x) && !is.matrix(y)) return(t(t(x) * y))
  if (is.matrix(y) && !is.matrix(x)) return(t(x * t(y)))
  stop('two arguments must contain exactly one matrix and one vector')
}
#' @rdname sweep_arithmetic
#' @export
`%_/_%` <- function(x, y) {
  if (is.matrix(x) && !is.matrix(y)) return(t(t(x) / y))
  if (is.matrix(y) && !is.matrix(x)) return(t(x / t(y)))
  stop('two arguments must contain exactly one matrix and one vector')
}

#' Check whether input is string
#' @param x object
#' @return boolean for whether x is string
#' @export
isString <- function(x){
  is.character(x) && (length(x) == 1)
}

#' assign colours to distinct values of a vector
#' @param v vector to be converted to colours
#' @return a vector of colours of equal length with two attributes: 'palette' giving the palette used for all distinct colours; 'uniq': vector of distinct values in v
#' @export
colorise <- function(v){
  uniq <- unique(v)
  ind <- match(v, uniq)
  palet <- matplotlib_palette(length(uniq))
  col <- palet[ind]
  attr(col, 'palette') <- palet
  attr(col, 'legend') <- uniq
  return(col)
}

#' assign line types to distinct values of a vector
#' @param v vector to be converted to line types / categories
#' @return a vector of line types of equal length with two attributes: 'lty' giving a vector of all line types used; 'uniq': vector of distinct values in v
#' @export
stylise <- function(v){
  uniq <- unique(v)
  ind <- match(v, uniq)
  lty <- seq_along(uniq)
  style <- lty[ind]
  attr(style, 'lty') <- lty
  attr(style, 'legend') <- uniq
  return(style)
}



#' generate a plot from a data frame with color and line type set by specific columns
#' @param x column name for the x values, needs to be a string
#' @param y column name for the y values, needs to be a string
#' @param col column name for the colour attributes, can be of any type, distinct values are represented by distinct colours
#' @param style column name for the line type attributes, can be of any type, distinct values are represented by distinct line types
#' @param data data frame
#' @param legend.position where to place the legend
#' @param ... other plotting parameters
#' @export
myplot <- function(x, y, col=NULL, style=NULL, data=.GlobalEnv, legend.position='topright', ...){
  xval <- data[[x]]
  yval <- data[[y]]
  legend.col <- legend.lty <- legend.txt <- c()

  if (is.null(col)) {
    plot.col <- legend.col <- 'black'
  } else {
    plot.col <- colorise(data[[col]])
    legend.col <- attr(plot.col, 'palet')
    legend.txt <- paste(col, '=', attr(plot.col, 'legend'))
  }
  if (is.null(style)){
    plot.lty <- legend.lty <- 1
  } else {
    plot.lty <- stylise(data[[style]])
    legend.lty <- attr(plot.lty, 'lty')
    legend.txt <- c(legend.txt, paste(style, '=', attr(plot.lty, 'legend')))
  }
  if (!is.null(col) && !is.null(style)){
    tmp1 <- legend.col; tmp2 <- legend.lty
    legend.col <- c(tmp1, rep_along('black', tmp2))
    legend.lty <- c(rep_along(1, tmp1), tmp2)
  }
  tmp <- data.frame(xval=xval, yval=yval, col=plot.col, lty=plot.lty)
  agg <- aggregate(cbind(xval, yval) ~ col + lty, data=tmp, FUN=list)
  plot(range(xval), range(yval), xlab=x, ylab=y, type='n', ...)
  for (i in 1:nrow(agg)){
    points(agg[i, 3][[1]], agg[i, 4][[1]], col=agg[i, 1], lty=agg[i, 2], type='b', ...)
  }
  legend(legend.position, legend=legend.txt, col=legend.col, lty=legend.lty)
}
