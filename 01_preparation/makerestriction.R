library(readr)
library(dplyr)
library(stringr)
library(tidyr)

# setwd("./000_data_IO/")

#読み込み：コーディング表
# fin <- "ch_coding_cellshape_motility_sporulation_ALL.txt"
fin_full <- "ch_coding_GTDB207_cellshape_motility_sporulation_ALLphylum.tsv"
mydata <- read_tsv(fin_full)

#ファイル名取得：門ごとの特徴コーディング表
# mykey <- "^ch_cellshape_motility_sporulation_"
mykey <- "^ch_GTDB207_cellshape_motility_sporulation_"
fin <- list.files(getwd(), pattern=mykey, full.names=F)

# #形質の抽出（ファイル名から）
# buf <- fin_full %>% strsplit(".tsv") %>% unlist
# buf <- buf %>% strsplit("_") %>% unlist
# #ALLを除去
# buf <- buf[-length(buf)]
# #GTDB207を除去
# buf <- buf[-1]
# mytraits <- buf[3:length(buf)]

#形質情報
mytraits <- c("cellshape", "motility", "sporulation")

#コーディング表の形質状態ごとに分割する
mybuf <- mydata %>% tidyr::separate(trait, into=mytraits, sep="_") %>% select(-ch)
mydata <- bind_cols(mydata, mybuf)

#門ごとに
#門ごとに存在する形質状態が異なり、登場するqAB = qZXなどが異なる
for (cfin in fin){
  #門名の抽出
  mystart <- regexpr(mykey, cfin) + nchar(mykey) - 1
  myend <- regexpr(".tsv", cfin) - 1
  mytaxa <- str_sub(cfin, mystart, myend)
  
  #読み込み：特徴コーディング表
  mydf_taxa <- read.table(cfin, header=F)
  colnames(mydf_taxa) <- c("accession", "ch")
  #形質状態を抽出
  mych <- mydf_taxa %>% distinct(ch) %>% arrange(ch) %>% pull(ch)
  
  #ファイル名
  fout <- fin_full %>% str_replace("ch", "mycommand") %>% str_replace("ALLphylum", mytaxa)
  
  #形質状態が２以上のもの
  if (length(mych) >= 2){
    #全部の組み合わせの文字生成
    mycomb <- combn(mych, 2) %>% as.matrix(ncol=2) %>% t %>% as.data.frame
    colnames(mycomb) <- c("ch1", "ch2")
    #列を逆さにしたものも生成し、加える
    buf <- mycomb %>% select(ch2, ch1)
    colnames(buf) <- c("ch1", "ch2")
    mycomb <- bind_rows(mycomb, buf)
    
    #形質の組み合わせに、形質ごとの列をマージ
    aaa <- left_join(mycomb, mydata, by=c("ch1" = "ch"))
    colnames(aaa)[3:ncol(aaa)] <- str_c("state1_", colnames(aaa)[3:ncol(aaa)], sep="")
    mydf <- left_join(aaa, mydata, by=c("ch2" = "ch"))
    colnames(mydf)[(ncol(aaa)+1):ncol(mydf)] <- str_c("state2_", colnames(mydf)[(ncol(aaa)+1):ncol(mydf)], sep="")
    
    #rate parametersの初期値（すべてに0を入力）
    #cellshape_motility_sporulationのうち２形質が同じものを選抜しmytag_rateをつけていく
    #そうでないものは、２形質が同時に変化しているので、初期値のrate=0でよい
    mydf$myrow <- seq(1, nrow(mydf))
    mydf$rate <- 0
    
    #形質の状態
    mys_cellshape <- mydf %>% distinct(state1_cellshape) %>% pull
    mys_motility <- mydf %>% distinct(state1_motility) %>% pull
    mys_sporulation <- mydf %>% distinct(state1_sporulation) %>% pull
    
    #Sporulationの遷移率に識別番号を入れる
    #言い換えれば、cellshapeは同じ、motilityも同じ行を抽出する
    mytag_rate <- 1
    for (cspo in mys_sporulation){
      for (ccs in mys_cellshape){
        #cell_shapeが同じもの選抜
        buf1 <- mydf %>% filter(state1_cellshape == ccs) %>% filter(state2_cellshape == ccs)
        #buf1がある場合
        if (nrow(buf1) > 0){
          #motilityが同じもの選抜
          for (cm in mys_motility){
            buf2 <- buf1 %>% filter(state1_motility == cm) %>% filter(state2_motility == cm)
            #buf2がある場合
            if (nrow(buf2) > 0){
              aaa <- buf2 %>% filter(state1_sporulation == cspo) %>% select(myrow) %>% pull
              mydf[aaa,]$rate <- mytag_rate
            }
          }
        }
      }
      mytag_rate <- mytag_rate + 1
    }
    
    #motilityの遷移率に識別番号を入れる
    for (cm in mys_motility){
      for (ccs in mys_cellshape){
        #cell_shapeが同じもの選抜
        buf1 <- mydf %>% filter(state1_cellshape == ccs) %>% filter(state2_cellshape == ccs)
        #buf1がある場合
        if (nrow(buf1) > 0){
          #sporulatioが同じもの選抜
          for (cspo in mys_sporulation){
            buf2 <- buf1 %>% filter(state1_sporulation == cspo) %>% filter(state2_sporulation == cspo)
            #buf2がある場合
            if (nrow(buf2) > 0){
              aaa <- buf2 %>% filter(state1_motility == cm) %>% select(myrow) %>% pull
              mydf[aaa,]$rate <- mytag_rate
            }
          }
        }
      }
      mytag_rate <- mytag_rate + 1
    }
    
    #cell_shapeの遷移率に識別番号を入れる
    for (ccs1 in mys_cellshape){
      for (ccs2 in mys_cellshape){
        for (cm in mys_motility){
          #cell_shapeが同じもの選抜
          buf1 <- mydf %>% filter(state1_motility == cm) %>% filter(state2_motility == cm)
          #buf1がある場合
          if (nrow(buf1) > 0){
            #sporulatioが同じもの選抜
            for (cspo in mys_sporulation){
              buf2 <- buf1 %>% filter(state1_sporulation == cspo) %>% filter(state2_sporulation == cspo)
              #buf2がある場合
              if (nrow(buf2) > 0){
                aaa <- buf2 %>% filter(state1_cellshape == ccs1) %>% filter(state2_cellshape == ccs2) %>%select(myrow) %>% pull
                #aaaがある場合
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
    
    #qABの生成
    mydf <- mydf %>% mutate(ch = str_c("q", ch1, ch2, sep=""))
    
    #rateが同じ数値をもつqAB,qBCをRestrictionコマンドで書き出す
    #rate = 0のとき
    aaa <- mydf %>% filter(rate == 0) %>% select(ch) %>% pull
    if (length(aaa) > 0){
      aaa <- aaa %>% paste(collapse=" ")
      mycommand <- str_c("Restrict", aaa, 0, sep=" ")
      #保存
      write.table(mycommand, fout, quote=F, append=T, col.names=F, row.names=F, sep="\t")
    }
          
    #それ以外のrateのとき
    #rateが同じ数値をもつqAB,qBCが２個以上存在している必要がある
    myratelist <- mydf %>% filter(rate != 0) %>% distinct(rate) %>% pull
    for (cr in myratelist){
      aaa <- mydf %>% filter(rate == cr) %>% select(ch) %>% pull
      #rateが２つ以上ある場合のみ
      if (length(aaa) > 1){
        aaa <- aaa %>% paste(collapse=" ")
        mycommand <- str_c("Restrict", aaa, sep=" ")
        write.table(mycommand, fout, quote=F, append=T, col.names=F, row.names=F, sep="\t")
      }
    }
  }
}

