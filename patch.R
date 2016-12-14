JSTOR_unpack1grams<-function (parallel = FALSE, path = getwd()) 
{
  if (substr(path, nchar(path), nchar(path)) == "/") {
    path <- substr(path, 1, nchar(path) - 1)
  }
  else {
    path <- path
  }
  setwd(paste0(path, "/wordcounts"))
  message("reading 1-grams into R...")
  myfiles <- dir(pattern = "\\.(csv|CSV)$", full.names = TRUE)
  suppressMessages(library(data.table))
  library(plyr)
  read_csv2dt <- function(x) data.table(fread(x, sep = ",", 
                                              stringsAsFactors = FALSE))
  if (parallel) {
    suppressMessages(library(snow))
    suppressMessages(library(parallel))
    cl <- makeCluster(detectCores(), type = "SOCK")
    clusterExport(cl, c("myfiles", "read_csv2dt"), envir = environment())
    clusterEvalQ(cl, library(data.table))
    aawc <- parLapplyLB(cl, myfiles, read_csv2dt)
    stopCluster(cl)
    invisible(gc(verbose = FALSE))
  }
  else {
    library(plyr)
    aawc <- llply(myfiles, read_csv2dt, .progress = "text", 
                  .inform = FALSE)
  }
  names(aawc) <- myfiles
  message("done")
  lens <- sapply(aawc, function(i) i[1]$WEIGHT + i[2]$WEIGHT + 
                   i[3]$WEIGHT)
  full <- unname(!is.na(lens))
  aawc1 <- aawc[full]
  message("reshaping the 1-grams into a document term matrix...")
  library(slam)
  library(tm)
  my_dtm_1gram <- function(x) {
    y <- as.integer(x$WEIGHT)
    names(y) <- x$WORDCOUNTS
    v = unname(y)
    i = rep(1, length(y))
    j = seq(1:length(y))
    z <- simple_triplet_matrix(v = v, i = i, j = j, nrow = max(i), 
                               ncol = max(j), dimnames = list(Docs = deparse(substitute(x)), 
                                                              Terms = names(y)))
    zz <- as.DocumentTermMatrix(z, weighting = weightTf)
    return(zz)
  }
  if (parallel) {
    suppressMessages(library(snow))
    suppressMessages(library(parallel))
    cl <- makeCluster(detectCores(), type = "SOCK")
    clusterExport(cl, c("aawc", "my_dtm_1gram"), envir = environment())
    clusterEvalQ(cl, list(library(tm), library(slam)))
    aawc1 <- parLapplyLB(cl, 1:length(aawc1), function(i) my_dtm_1gram(aawc1[[i]]))
    stopCluster(cl)
    rm("aawc")
    invisible(gc(verbose = FALSE))
  }
  else {
    library(plyr)
    aawc2 <- llply(1:length(aawc1), function(i) my_dtm_1gram(aawc1[[i]]), 
                   .progress = "text", .inform = FALSE)
  }
  myfiles1 <- myfiles[full]
  library(stringr)
  names(aawc2) <- str_extract(basename(myfiles1), "[^wordcounts_].+[^.CSV]")
  message("done")
  message("arranging bibliographic data...")
  setwd(path)
  read_citations <- function(i) {
    if (stringr::str_sub(i, start = -3) == "CSV" | stringr::str_sub(i, 
                                                                    start = -3) == "csv") {
      read.csv(i, quote = "", row.names = NULL, comment.char = "", 
               header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
    }
    else {
      if (stringr::str_sub(i, start = -3) == "TSV" | stringr::str_sub(i, 
                                                                      start = -3) == "tsv") {
        read.delim(i, row.names = NULL, comment.char = "", 
                   header = TRUE, stringsAsFactors = FALSE, colClasses = "character", 
                   quote = "")
      }
      else {
        "Citations files cannot be loaded"
      }
    }
  }
  cit <- read_citations(dir(pattern='citations'))
  library(stringr)
  cit$id <- str_extract(chartr("/", "_", cit$id), ".*[^\t]")
  citfla <- cit[cit$publisher == "fla", ]
  citfla$id <- as.character(gsub(" ", "", citfla$id))
  wordcounts <- aawc2[which(names(aawc2) %in% citfla$id)]
  bibliodata <- (merge(names(wordcounts), citfla, by.x = 1, 
                       by.y = "id"))
  bibliodata$year <- str_extract(bibliodata$issue, "[[:digit:]]{4}")
  rm(aawc1, aawc2, cit, citfla, myfiles)
  invisible(gc(verbose = FALSE))
  if (parallel) {
    suppressMessages(library(snow))
    suppressMessages(library(parallel))
    cl <- makeCluster(detectCores(), type = "SOCK")
    clusterExport(cl, c("wordcounts"), envir = environment())
    clusterEvalQ(cl, library(tm))
    wordcounts <- do.call(tm:::c.DocumentTermMatrix, wordcounts)
    stopCluster(cl)
    invisible(gc(verbose = FALSE))
  }
  else {
    wordcounts <- do.call(tm:::c.DocumentTermMatrix, wordcounts)
  }
  wordcounts$dimnames$Docs <- as.character(bibliodata$x)
  wordcounts <- wordcounts[unique(as.character(wordcounts$dimnames$Docs[1:nrow(wordcounts)])), 
                           ]
  message("removing stopwords...")
  wordcounts <- wordcounts[, !(wordcounts$dimnames$Terms %in% 
                                 stopwords(kind = "en"))]
  message("done")
  message("discarding words with <3 characters (probably OCR errors)...")
  wordcounts <- wordcounts[, nchar(wordcounts$dimnames$Terms) > 
                             3]
  message("done")
  message("discarding words with >2 consecutive characters (probably OCR errors)...")
  wordcounts <- wordcounts[, !grepl("(.)\\1{2,}", wordcounts$dimnames$Terms)]
  message("done")
  message("discarding non-ASCII characters...")
  wordcounts <- wordcounts[, (wordcounts$dimnames$Terms %in% 
                                iconv(wordcounts$dimnames$Terms, "latin1", "ASCII", sub = ""))]
  message("done")
  message("finished with 1-grams")
  return(list(wordcounts = wordcounts, bibliodata = bibliodata))
}
