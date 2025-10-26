# このスクリプトは、mycommand(BayesTraitsの制約を書いたファイル)に書かれた遷移だけを考慮して、Model stringのどの列がどの遷移率に相当するかを判断している。
# そのため、例えば、制約に書かれていない遷移（要素数が多いときには、同じ要素の獲得・消失が１遷移だけの場合がある）は考慮していない。

library(readr)
library(dplyr)
library(stringr)
library(tidyr)

#パスの読み込み
#ファイル名のみ（フルパスではない）
# setwd("./000_data_IO")
# setwd("../000_data_IO_butterfly/")
fin <- list.files(getwd(), pattern="^eData_", full.names=F)
#independent解析を削除
buf <- grep("_i_", fin)
if (length(buf) > 0){
  fin <- fin[-buf]
}

#読み込み：制限を加えたパラメーター
fin_res <- list.files(getwd(), pattern="^mycommand_coding_", full.names=F)

for (cfin in fin){
  print(cfin)
  #形質＋分類群情報の抽出
  #mycommandファイルとの紐づけに利用する
  buf <- str_replace(cfin, "eData_ch_", "")
  mykey <- str_replace(buf, ".tsv.tsv", "")
  
  #読み込み
  mydata <- read_tsv(cfin)

  #Model string列を抽出
  mycol <- grep("string", colnames(mydata))
  buf <- mydata %>% pull(mycol)
  #Model string列の文字数をカウント
  mynum_string <- buf[1] %>% as.character %>% nchar
  #１文字目の'を除去して抽出
  mydf <- str_sub(pull(mydata, mycol), start = 2, end = mynum_string) %>% as.matrix(ncol=1) %>% as.data.frame()
  colnames(mydf) <- "datacol"
  
  #遷移率の列名を取得：qAB,qAC,...
  mystart <- grep("^q", colnames(mydata))[1]
  myend <- grep("^Root", colnames(mydata))[1] - 1

  myrates <- c(colnames(mydata)[mystart:myend])
  
  #遷移率の制限がある場合
  #
  if (length(fin_res) > 0){
    #遷移率の制限の読み込み
    mykey2 <- str_c(mykey, ".tsv", sep="")
    fin_res_taxon <- grep(mykey2, fin_res)
    if (length(fin_res_taxon) > 0){
      mydata_res <- readLines(fin_res[fin_res_taxon])
      
      #遷移率＝0にした遷移率
      buf <- grep("0", mydata_res)
      bufbuf <- mydata_res[buf]
      bufbufbuf <- str_split(bufbuf, " ") %>% unlist
      buf4 <- grep("q", bufbufbuf)
      myres_zero <- bufbufbuf[buf4]
      
      #同じ値をとる遷移率の内、最初以外を集める
      buf <- grep("0", mydata_res)
      mydata_res_mod <- mydata_res[-buf]
      myres_all <- myres_zero
      for (i in 1:length(mydata_res_mod)){
        aaa <- str_split(mydata_res_mod[i], " ") %>% unlist
        buf <- grep("q", aaa)
        bufbuf <- aaa[buf] %>% sort
        myres_all <- c(myres_all, bufbuf[-1])
      }
      
      #以下は何のためにやってるの？
      #遷移率を固定するとModel stringに出てこないが、qABのところには固定した遷移率が表示される。
      #そこで、eDataの列名リストから、制限して=0にしたものを削除する。
      myrates <- setdiff(myrates, myres_all)
      
      #以下は記録のために残しているだけ
      #以下のやり方はダメ。制約してなくても=0になってしまう場合も削除してしまう。
      # #そこでqABの値をuniqueして、固定した遷移率＝nrow1なら削除している。
      # for (i in length(myrates):1){
      #   buf <- mydata[myrates[i]] %>% unique %>% nrow
      #   if (buf == 1){
      #     myrates <- myrates[-i]
      #   }
      # }
    }
  }
  
  myinto <- myrates
  
  #Model string列を１文字１列にする
  #スペース区切りでsplit、新しい列にする
  mydf_z <- mydf %>% tidyr::separate(datacol, into = myinto, sep = " ")

  #各列ごと(遷移率ごとqAB,qAC,...)に"Z"の個数を数える
  mycount_z <- NULL
  for (i in 1:ncol(mydf_z)){
    mycount_z[i] <- mydf_z[mydf_z[,i]=="Z",i] %>% length
  }
  
  #百分率に変換
  mypercent <- mycount_z * 100 / nrow(mydf_z)
  
  #出力へデータ整形
  #１列目：遷移、２列目：z-score
  #１列目：遷移の作成
  myinto <- myinto %>% as.matrix(ncol = 1) %>% as.data.frame
  colnames(myinto) <- "parameter"
  #２列目：z-scoreの作成
  mypercent <- mypercent %>% as.matrix(ncol = 1) %>% as.data.frame
  colnames(mypercent) <- "zscore"
  #１列目、２列目の統合
  myresult <- bind_cols(myinto, mypercent)
  
  #遷移率の制限がある場合
  if (length(fin_res) > 0 && length(fin_res_taxon) > 0){
    #遷移率＝０の遷移率を追加する
    if (length(myres_zero) > 0){
      mydf <- data.frame(parameter=myres_zero, zscore=100)
      myresult <- bind_rows(myresult, mydf)
    }
    
    #同じ値をとる遷移率を追加する
    buf <- grep("0", mydata_res)
    mydata_res_mod <- mydata_res[-buf]
    if (length(mydata_res_mod) > 0){
      for (i in 1:length(mydata_res_mod)){
        aaa <- str_split(mydata_res_mod[i], " ") %>% unlist
        buf <- grep("q", aaa)
        bufbuf <- aaa[buf] %>% sort
        myvalue <- myresult %>% filter(parameter == bufbuf[1]) %>% pull(zscore)
        mydf <- data.frame(parameter=c(bufbuf[-1]), zscore=myvalue)
        myresult <- bind_rows(myresult, mydf)
      }
    }
  }
  
  #ファイル名
  fout <- str_replace(cfin, "eData", "zscore")
  write_tsv(myresult, fout)
}

