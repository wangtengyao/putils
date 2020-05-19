########## string operations ##########

#' String length
#' @param str a string
#' @return its length
#' @export
strlen = function(str){
  return(nchar(str))
}

#' Substring search
#' @description find the position of the nth occurrence of needle in haystack, returns 0 if not found
#' @param haystack a string
#' @param needle substring to search for
#' @param startpos start position for search
#' @param n the nth occurrence
#' @return an integer
#' @export
strstr <- function(haystack, needle, startpos=1, n=1){
  aa <- unlist(strsplit(substring(haystack, startpos), needle))
  if(length(aa) < n + 1 ) return(0);
  return(sum(nchar(aa[1:n])) + startpos + (n-1) * nchar(needle) )
}

#' Extract left substring
#' @param str a string
#' @param len length of substring
#' @export
#'
prefix <- function(str, len){
  substr(str, 1, len)
}

#' Extract right substring
#' @param str a string
#' @param len length of substring
#' @export
#'
suffix <- function(str, len){
  substr(str, strlen(str) - len + 1, strlen(str))
}

