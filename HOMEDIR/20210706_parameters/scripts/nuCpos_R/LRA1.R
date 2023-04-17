LRA1 <- function(inseq, species = "mm", silent = FALSE){

    if(silent == FALSE) message("species: ", species, "\n")
    if(!is(inseq)[1] == "character"){
        if(is(inseq)[1] == "DNAString"){
            if(requireNamespace("Biostrings", quietly = TRUE)){
                inseq <- as.character(inseq)
                if(silent == FALSE){
                message(
                "The class of inseq was changed from DNAString to character")
                }
            }else{
            message("DNAString cannot be changed to a character string")
            out <- NA; names(out) <- "LRA"; return(out)
            }
        }else{
            message("The class of inseq must be DNAString or character")
            out <- NA; names(out) <- "LRA"; return(out)
        }
    }
    if(nchar(inseq) != 147){
        message("Length of inseq: ", nchar(inseq), "bp")
        message("The length of inseq must be 147 bp")
        out <- NA; names(out) <- "LRA"; return(out)
    }

    if(species == "sc"){
        freqL1 <- nature11142_s2.147.freqL
        tranL1 <- nature11142_s2.147.tranL
        TtranL2 <- nature11142_s2.147.tranL2
        TtranL3 <- nature11142_s2.147.tranL3
        TtranL4 <- nature11142_s2.147.tranL4
        TfreqN4L <- nature11142_s2.147.freqN4SA
        TfreqN4R <- nature11142_s2.147.freqN4R
        TtranN4 <- nature11142_s2.147.tranN4_SMA
    }

    if(species == "scH3H4_33L"){
        freqL1 <- nature11142_s2.147.freqL
        tranL1 <- nature11142_s2.147.tranL
        TtranL2 <- nature11142_s2.147.tranL2
        TtranL3 <- nature11142_s2.147.tranL3
        TtranL4 <- nature11142_s2.147.tranL4
        TfreqN4L <- merged33_67352.147.freqN4SA
        TfreqN4R <- merged33_67352.147.freqN4R
        TtranN4 <- merged33_67352.147.tranN4_SMA
    }

    inseq_num <- as.integer(charToRaw(inseq))
    outlist <- .Fortran("LRA_1", 
            inseq_num = inseq_num,
            logascL = numeric(length=1), logascR = numeric(length=1), 
            freqL1 = freqL1, 
            tranL1 = tranL1, TtranL2 = TtranL2, TtranL3 = TtranL3, 
            TtranL4 = TtranL4, TfreqN4L = TfreqN4L, TfreqN4R = TfreqN4R, 
            # TtranN4 = TtranN4, PACKAGE = "nuCpos")[c(2,3)]
            TtranN4 = TtranN4)[c(2,3)]

    # N を含んでいるときはすべてを NA にする。
    if(outlist[[1]] == 0){
        outlist[1] <- as.numeric(NA)    
    }
    out <- unlist(outlist)
    names(out) <- c("L", "R")
    return(out)
}


