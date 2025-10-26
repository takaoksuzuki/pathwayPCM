library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(igraph)
library(RColorBrewer)
library(Polychrome)

# setwd("./000_data_IO/")

#読み込み：z-scoreの結果のファイル名
#-----------------------------------------------------------------------------------------
fin_z <- list.files(getwd(), pattern ="^zscore_ch_GTDB207_cellshape_motility_sporulation_", full.names = F)
#-----------------------------------------------------------------------------------------

#読み込み：qrateの結果
#for-loop内で読み込んでいる

for (cfin in fin_z) {
  #ファイル名の生成
  fin_q <- str_replace(cfin, "zscore", "qrate")

  #読み込み：z-score
  #read_tsvにすると色を読み込まないので注意！！！
  mydata_z <- read.delim(cfin)
  #列名
  colnames(mydata_z) <- c("parameter","zscore")

  #読み込み：qrates+GLMでの分類
  #read_tsvにすると色を読み込まないので注意！！！
  mydata_q <- read.delim(fin_q)
  #列名
  colnames(mydata_q) <- c("parameter","qrates")
  
  #マージ
  mydata <- left_join(mydata_z, mydata_q, by="parameter")
  
  #エッジリストの作成
  #状態の分割
  #qAB -> A, B
  mybuf <- mydata %>% mutate(state_1=str_sub(parameter, 2, 2)) 
  mydf <- mybuf %>% mutate(state_2=str_sub(parameter, 3, 3)) %>% select(state_1, state_2, zscore, qrates)
  
  #変換：z-score -> non-z-score
  mydf <- mydf %>% mutate(nonzscore = 100 - zscore) %>% select(-zscore)
  
  #閾値
  #-------------------------------
  mythreshold_lower <- 30
  #-------------------------------
  #閾値で切る
  myedge <- mydf %>% filter(nonzscore > mythreshold_lower)

  if (nrow(mydf) > 0){
    #グラフに変換
    g <- graph_from_data_frame(myedge, directed = TRUE)
    # #重み
    E(g)$weight <- myedge$qrates

    #描画
    #----------------------------------------------------------
    #描画レイアウト
    mylayout=layout_with_graphopt

    #ファイル名
    fout <- str_replace(cfin, "zscore_ch", "plot")
    fout <- str_replace(fout, "tsv", "png")
    
    #保存
    png(fout, width=3000, height=2000)
    
    par(mar = c(0, 0, 0, 0)) #  余白の広さを行数で指定．下，左，上，右の順．
    # ---- プロット ----
    plot.igraph(g,
                layout = mylayout,
                vertex.color = "white",
                vertex.size = 10,
                vertex.shape = "circle",
                vertex.label.cex = 5,
                vertex.frame.width=5,
                edge.lty = E(g)$lty,
                edge.arrow.size = 4,
                edge.curved = 0.2)
    dev.off()
  }
}
