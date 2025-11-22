library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# === Set data directory ===
data_dir <- file.path("data")

# === Load coding table ===
fin_full <- file.path(data_dir, "ch_coding_GTDB207_cellshape_motility_sporulation_ALLphylum.tsv")
mydata <- readr::read_tsv(fin_full)

# === Get file names: trait coding tables by phylum ===
mykey <- "^ch_GTDB207_cellshape_motility_sporulation_"
fin <- list.files(data_dir, pattern = mykey, full.names = TRUE)

# === Define target traits ===
mytraits <- c("cellshape", "motility", "sporulation")

# === Split the coding table by trait states ===
mybuf <- mydata %>% tidyr::separate(trait, into=mytraits, sep="_") %>% select(-ch)
mydata <- bind_cols(mydata, mybuf)

# === Loop over each phylum ===
# Each phylum may have different sets of trait states, resulting in distinct qAB (e.g., qZX) transitions.
for (cfin in fin){
  # Extract phylum name from file name
  mystart <- regexpr(mykey, cfin) + nchar(mykey) - 1
  myend <- regexpr(".tsv", cfin) - 1
  mytaxa <- str_sub(cfin, mystart, myend)
  
  # === Read the trait coding table ===
  mydf_taxa <- read.table(cfin, header=F)
  colnames(mydf_taxa) <- c("accession", "ch")
  # Extract trait states
  mych <- mydf_taxa %>% distinct(ch) %>% arrange(ch) %>% pull(ch)
  
  # output file name
  fout <- fin_full %>% str_replace("ch", "mycommand") %>% str_replace("ALLphylum", mytaxa)
  
  # === Process only if two or more trait states exist ===
  if (length(mych) >= 2){
    # Generate all possible pairwise combinations of trait states
    mycomb <- combn(mych, 2) %>% as.matrix(ncol=2) %>% t %>% as.data.frame
    colnames(mycomb) <- c("ch1", "ch2")

    # Add reversed combinations (e.g., AB and BA)
    buf <- mycomb %>% select(ch2, ch1)
    colnames(buf) <- c("ch1", "ch2")
    mycomb <- bind_rows(mycomb, buf)
    
    # Merge each trait combination with the corresponding trait columns
    aaa <- left_join(mycomb, mydata, by=c("ch1" = "ch"))
    colnames(aaa)[3:ncol(aaa)] <- str_c("state1_", colnames(aaa)[3:ncol(aaa)], sep="")
    mydf <- left_join(aaa, mydata, by=c("ch2" = "ch"))
    colnames(mydf)[(ncol(aaa)+1):ncol(mydf)] <- str_c("state2_", colnames(mydf)[(ncol(aaa)+1):ncol(mydf)], sep="")
    
    # === Initialize rate parameters (set all to 0) ===
    # For each of the three traits (cellshape, motility, sporulation),
    # assign rate IDs (mytag_rate) to transitions where only one trait changes.
    # Transitions involving simultaneous changes in two traits remain rate = 0.
    mydf$myrow <- seq(1, nrow(mydf))
    mydf$rate <- 0
    
    # === Extract possible trait states ===
    mys_cellshape <- mydf %>% distinct(state1_cellshape) %>% pull
    mys_motility <- mydf %>% distinct(state1_motility) %>% pull
    mys_sporulation <- mydf %>% distinct(state1_sporulation) %>% pull
    
    # === Assign rate IDs for sporulation transitions ===
    # In other words, select rows where both cellshape and motility remain constant.
    mytag_rate <- 1
    for (cspo in mys_sporulation){
      for (ccs in mys_cellshape){
        # Select rows where cellshape is identical
        buf1 <- mydf %>% filter(state1_cellshape == ccs) %>% filter(state2_cellshape == ccs)
        # If such rows exist
        if (nrow(buf1) > 0){
          # Select rows where motility is identical
          for (cm in mys_motility){
            buf2 <- buf1 %>% filter(state1_motility == cm) %>% filter(state2_motility == cm)
            # If such rows exist
            if (nrow(buf2) > 0){
              aaa <- buf2 %>% filter(state1_sporulation == cspo) %>% select(myrow) %>% pull
              mydf[aaa,]$rate <- mytag_rate
            }
          }
        }
      }
      mytag_rate <- mytag_rate + 1
    }
    
    # === Assign rate IDs for motility transitions ===
    for (cm in mys_motility){
      for (ccs in mys_cellshape){
        # Select rows where cellshape is identical
        buf1 <- mydf %>% filter(state1_cellshape == ccs) %>% filter(state2_cellshape == ccs)
        # If such rows exist
        if (nrow(buf1) > 0){
          # Select rows where sporulation is identical
          for (cspo in mys_sporulation){
            buf2 <- buf1 %>% filter(state1_sporulation == cspo) %>% filter(state2_sporulation == cspo)
            # If such rows exist
            if (nrow(buf2) > 0){
              aaa <- buf2 %>% filter(state1_motility == cm) %>% select(myrow) %>% pull
              mydf[aaa,]$rate <- mytag_rate
            }
          }
        }
      }
      mytag_rate <- mytag_rate + 1
    }
    
    # === Assign rate IDs for cellshape transitions ===
    for (ccs1 in mys_cellshape){
      for (ccs2 in mys_cellshape){
        for (cm in mys_motility){
          # Select rows where motility is identical
          buf1 <- mydf %>% filter(state1_motility == cm) %>% filter(state2_motility == cm)
          # If such rows exist
          if (nrow(buf1) > 0){
            # Select rows where sporulation is identical
            for (cspo in mys_sporulation){
              buf2 <- buf1 %>% filter(state1_sporulation == cspo) %>% filter(state2_sporulation == cspo)
              # If such rows exist
              if (nrow(buf2) > 0){
                aaa <- buf2 %>% filter(state1_cellshape == ccs1) %>% filter(state2_cellshape == ccs2) %>%select(myrow) %>% pull
                # If matches found
                if (length(aaa) > 0){
                  mydf[aaa,]$rate <- mytag_rate
                }
              }
            }
          }
        }
        mytag_rate <- mytag_rate + 1
      }
    }
    
    # === Generate qAB labels ===
    mydf <- mydf %>% mutate(ch = str_c("q", ch1, ch2, sep=""))
    
    # === Write Restrict commands for rate groups ===
    # Case 1: rate = 0
    aaa <- mydf %>% filter(rate == 0) %>% select(ch) %>% pull
    if (length(aaa) > 0){
      aaa <- aaa %>% paste(collapse=" ")
      mycommand <- str_c("Restrict", aaa, 0, sep=" ")
      # Save to file
      write.table(mycommand, fout, quote=F, append=T, col.names=F, row.names=F, sep="\t")
    }
          
    # Case 2: nonzero rate groups
    # Only include groups that have two or more qAB transitions sharing the same rate ID
    myratelist <- mydf %>% filter(rate != 0) %>% distinct(rate) %>% pull
    for (cr in myratelist){
      aaa <- mydf %>% filter(rate == cr) %>% select(ch) %>% pull
      # Process only if two or more transitions exist
      if (length(aaa) > 1){
        aaa <- aaa %>% paste(collapse=" ")
        mycommand <- str_c("Restrict", aaa, sep=" ")
        write.table(mycommand, fout, quote=F, append=T, col.names=F, row.names=F, sep="\t")
      }
    }
  }
}



