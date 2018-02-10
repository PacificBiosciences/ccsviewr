library(grid)
mkDifString <- function(ref, ccs) {
  r2 = strsplit(ref, "")[[1]]
  c2 = strsplit(ccs, "")[[1]]
  matches = rep(".", length(c2))
  matches[r2!=c2] = "*"
  paste(matches, sep="", collapse = "")
}


#' Print a PDF of the MSA in a dataframe returned by the AlnsToDataFrame function.
#'
#'
#' @param df A data frame returned by a call to AlnsToDataFrame.
#' @param pdfname The name of the pdf file to output
#' @param start Optionally subset a portion of the alignment by specifying a start column.
#' @param end Optionally subset a portion of the alignment by specifying an end column.
#' @param useOrigColors Whether to use my original color scheme http://colorbrewer2.org/?type=qualitative&scheme=Set1&n=4
#'
#' @return Returns a grid object with the alignment
#' @export
plotMSA <- function(df, pdfname, start=0, end=-1, showPositions=TRUE, useOrigColors = FALSE, showBases=TRUE, showDifString = TRUE) {

  # Prepare and subset the data
  ids = as.character(df$id)
  seqs = as.character(df$seq)
  if (end < 0) {
    end = nchar(seqs[1])
  }
  seqs = sapply(seqs, function(x) substr(x, start, end))
  names(seqs) <- NULL

  # Add in difference string
  if(showDifString) {
    dif = mkDifString(seqs[1], seqs[2])
    seqs = c(dif, seqs)
    ids = c("", ids)
  }


  n_seq = nchar(seqs[1])
  if(showPositions) {
    step = 10
    ticks = seq(step, n_seq, step)
    res = lapply(ticks, function(x) paste(paste(rep("-", step - nchar(x) -1), sep="", collapse=""), "|", x, sep="", collapse="" ))
    seqs = c(paste(res, collapse=""), seqs)
    ids = c("", ids)
  }
  idlengths = sapply(ids, nchar)
  maxid = max(idlengths)
  longid = ids[which(idlengths==maxid)][1]
  maxseq = max(sapply(seqs, nchar))



  ## Now to setup the graphics. Note that this is very annoying because I need
  ## to know the absolute size in inches when calling pdf() but I won't know that
  ## until all the plotting is done. So instead will apply some heuristics to the
  ## size of strings as returned by the strwidth function

  # Set us up to output in courier
  old = par()
  on.exit(suppressWarnings(par(old)))
  par(family="Courier")
  yb = unit(1, "char")

  # Space factor is just the borders I want
  spaceFactor = .5
  # Hard measurements + Heuristics to get sizing in inches
  idWidth = strwidth(longid, units="inches")
  seqWidth = strwidth(seqs[1], units="inches") * 1.5
  N = length(seqs)
  totX = (idWidth + seqWidth) + spaceFactor
  totY = yb * N
  pdfY = as.numeric(convertUnit(totY, "inches")) + spaceFactor

  # Output PDF
  pdf(pdfname, width=totX, height = pdfY, family="Courier")

  # Print in grid layout
  grid.newpage()
  base = grid.layout(ncol=2, widths = unit(c(idWidth, seqWidth), units = "inches"))
  vpbase = viewport(layout=base)
  pushViewport(vpbase)
  grid.rect()
  ypositions = rev((1/(N)) *seq(0, N-1 ))
  xpositions = rev((1/maxseq) * seq(0, maxseq-1))
  just = c("left", "bottom")

  # Print the IDs
  vpids = viewport(layout.pos.col = 1, name="ids")
  pushViewport(vpids)
  idmat = grid.layout(nrow = N)
  vpidmat = viewport(layout=idmat)
  pushViewport(vpidmat)
  drawId <- function(i)
  {
    vpc = viewport(layout.pos.row = i)
    pushViewport(vpc)
    if(ids[i] == "dark") {
      grid.text(ids[i], just=c("right", "center"), x=unit(0.98, "npc") ,
                gp=gpar(fontface=2, col="black"))
    } else if (ids[i] == "light") {
      grid.text(ids[i], just=c("right", "center"), x=unit(0.98, "npc"),
                gp=gpar(fontface=3))#, y=ypositions[i], x=0, just=just)
    }
    else if (ids[i] == "bunnanda") {
      grid.text(ids[i], just=c("right", "center"), x=unit(0.98, "npc"),
                gp=gpar(col="#1f78b4", fontface=2))#, y=ypositions[i], x=0, just=just)
    }
    else {
      grid.text(ids[i], just=c("right", "center"), x=unit(0.98, "npc"), gp=gpar(col="#33a02c", fontface=2))
    }
    #377eb8
    popViewport()
  }
  lapply(1:N, drawId)
  popViewport()
  popViewport()

  # Print the aligment columns
  vpseqs = viewport(layout.pos.col = 2, name="seqs")
  pushViewport(vpseqs)
  mat = grid.layout(ncol = maxseq, nrow = N)#, respect = TRUE)
  vpmat = viewport(layout=mat)
  pushViewport(vpmat)
  orgFillColors = c("#33a02c", "#a6cee3", "#1f78b4", "#b2df8a")
  newFillColors = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
  newTextColors = rep("black", 4)
  orgTextColors = newTextColors
  orgTextColors[3] = "white"
  fillColors = newFillColors
  textColors = newTextColors
  if (useOrigColors) {
    fillColors = orgFillColors
    textColors = orgTextColors
  }
  drawSeq <- function(i) {
    seq = strsplit(seqs[i], "")[[1]]
    drawBP <- function(j) {
      vpc = viewport(layout.pos.col = j, layout.pos.row = i)
      pushViewport(vpc)
      # Color scheme from http://colorbrewer2.org/
      fcolor = "white"
      tcolor = "black"
      bp = seq[j]
      if(bp == "A") {
        fcolor = fillColors[1]
        tcolor = textColors[1]
      } else if(bp == "C") {
        fcolor = fillColors[2]
        tcolor = textColors[2]
      } else if (bp == "G") {
        fcolor = fillColors[3]
        tcolor = textColors[3]
      } else if(bp=="T") {
        fcolor = fillColors[4]
        tcolor = textColors[4]
      }
      grid.rect(gp = gpar(fill=fcolor, col=NA))
      if(showBases) {
        grid.text(bp, gp=gpar(col=tcolor, cex=0.9))
      }
      popViewport()
    }
    sapply(1:length(seq), drawBP)
  }
  sapply(1:N, drawSeq)
  popViewport()
  popViewport()
  dev.off()


}

#plotMSA(df, "temp.pdf")1
