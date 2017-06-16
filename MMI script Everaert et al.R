-+#install packages AICcmodavg
  #utils:::menuInstallLocal()
  require(AICcmodavg)
  library(AICcmodavg)
  require(verification)
  require(ROCR)
  require(pROC)
  require(mgcv) 

  head(dataframe)
  dim(dataframe)
  
  dataframe$Acartia <- dataframe$Acartia/max(dataframe$Acartia, na.rm = TRUE)
  
  dataframe$TEMP    <- dataframe$TEMP/max(dataframe$TEMP, na.rm = TRUE)
  dataframe$SAL     <- dataframe$SAL/max(dataframe$SAL, na.rm = TRUE)
  dataframe$SLCA    <- dataframe$SLCA/max(dataframe$SLCA, na.rm = TRUE)
  dataframe$DIN     <- dataframe$DIN/max(dataframe$DIN, na.rm = TRUE)
  dataframe$CHFLa   <- dataframe$CHFLa/max(dataframe$CHFLa, na.rm = TRUE)
  dataframe$PCB     <- dataframe$PCB/max(dataframe$PCB, na.rm = TRUE)
  dataframe$Diatomratio <- dataframe$Diatomratio/max(dataframe$Diatomratio, na.rm = TRUE)
  
  head(dataframe)
  dim(dataframe)

  n.predictors <- c("2","3","4","5")
  predictors <- c("s(CHFLa, k = 4)","s(TEMP, k = 4)","s(SAL, k = 4)","s(PCB, k = 4)", "s(Diatomratio, k = 4)")
  
  
  # Acartia one predictor  - make equations yourself - store models outputs, aic, aicc and drop these in a matrix
  response <- "log10(Acartia)"
  resultaten.Acartia <- NULL
  aics <- NULL
  aiccs <- NULL
  formules <- c("log10(Acartia) ~ s(CHFLa, k = 4)","log10(Acartia) ~ s(TEMP, k = 4)","log10(Acartia) ~ s(SAL, k = 4)","log10(Acartia) ~ s(PCB, k = 4)","log10(Acartia) ~ s(Diatomratio, k = 4)")
  models <- list()
  AICs   <- rep(NA,length(formules))
  AICcs   <- rep(NA,length(formules))
  all.roc  <- rep(NA,length(formules))
  
  coeff.all.1 <- NULL
  coeff.all.2 <- NULL
  coeff <- NULL
  
  for (i in 1:length(formules))
  {
    coeff <- NULL
    model <- NULL
    pred <- NULL
    roc <- NULL
    
    model  <- gam(as.formula(formules[i]),data=dataframe)
    models[[i]] <- model
    aic  <- NULL
    aicc <- NULL
    aic  <- AIC(model)
    aicc <- AICc(model)
    AICs[i]  <- aic
    AICcs[i] <- aicc
    
    # calculation of ROC
    pred<-predict(model,newdata=dataframe,type="response")
    datasetAcartia <- as.matrix(dataframe$Acartia)

    roc <- auc(as.numeric(datasetAcartia),as.numeric(pred))
    all.roc[i] <- roc
    
    coeff <- matrix(data=0,nrow=1,ncol=6)
    colnames(coeff) <- c("(Intercept)","s(CHFLa).3","s(TEMP).3","s(SAL).3","s(PCB).3","s(Diatomratio).3")
    
    matrix.coefficient <- as.matrix(model$coefficients)
    u <- 1
    
    for (u in 1:length(colnames(coeff)))
    {
      test <- NA
      test <- which(rownames(matrix.coefficient[])==colnames(coeff)[u])
      if (length(test)>0){
        coeff[1,colnames(coeff)[u]] <- matrix.coefficient[test,1] 
      }else{
        coeff[1,colnames(coeff)[u]] <- 0
      }
    }
    
    coeff.all.2 <- rbind(coeff.all.2,coeff)
  }
  resultaat.Acartia1 <- cbind(formules,models,AICs,AICcs,all.roc)
  resultaat.Acartia1 <- as.data.frame(resultaat.Acartia1)
  resultaat.Acartia1 <- as.matrix(resultaat.Acartia1)
  resultaat.Acartia1 <- resultaat.Acartia1[,-c(2,5)]
  
  
  # Acartia two or more predictors  - store models outputs, aic, aicc and drop these in a matrix
  for (p in 1:length(n.predictors))
  {
    coeff <- NULL
    models <- list()
    
    predictorCombs <- combn(predictors,as.numeric(n.predictors[p]))
    
    AICs   <- rep(NA,ncol(predictorCombs))
    AICcs   <- rep(NA,ncol(predictorCombs))
    all.roc   <- rep(NA,ncol(predictorCombs))
    
    forms <- paste(response," ~ ",predictorCombs[1,],sep="")
    
    for (o in c(2:as.numeric(n.predictors[p])))
    {
      forms <- paste(forms,"+",predictorCombs[o,],sep="")
    }
    
    for (i in c(1:length(forms)))
    {
      model <- NULL
      aic   <- NULL
      aicc   <- NULL
      roc <- NULL
      
      form <- forms[i]
      model  <- gam(as.formula(form),data=dataframe)
      models[[i]] <- model
      aic <- AIC(model)
      aicc <- AICc(model)
      AICs[i] <- aic
      AICcs[i] <- aicc
      
      pred<-predict(model,newdata=dataframe,type="response")
      datasetAcartia <- as.matrix(dataframe$Acartia)
      
      roc <- auc(as.numeric(datasetAcartia),as.numeric(pred))
      all.roc[i] <- roc
      
      
      coeff <- matrix(data=0,nrow=1,ncol=6)
      colnames(coeff) <- c("(Intercept)","s(CHFLa).3","s(TEMP).3","s(SAL).3","s(PCB).3","s(Diatomratio).3")
      matrix.coefficient <- as.matrix(model$coefficients)
      
      u<-1
      
      for (u in 1:length(colnames(coeff)))
      {
        test <- NA
        test <- which(rownames(matrix.coefficient[])==colnames(coeff)[u])
        if (length(test)>0){
          coeff[1,colnames(coeff)[u]] <- matrix.coefficient[test,1] 
        }else{
          coeff[1,colnames(coeff)[u]] <- 0
        }
      }
      coeff.all.1 <- rbind(coeff.all.1,coeff)
    }
    
    resultaat.Acartia<-cbind(forms,models,AICs,AICcs,all.roc)
    resultaten.Acartia <- rbind(resultaten.Acartia,resultaat.Acartia)
  }
  resultaten.Acartia <- as.data.frame(resultaten.Acartia)
  resultaten.Acartia <- as.matrix(resultaten.Acartia)
  resultaten.Acartia <- resultaten.Acartia[,-c(2,5)]
  
  resultaten.Acartia <- rbind(resultaat.Acartia1,resultaten.Acartia)
  
  #write.csv(resultaten.Acartia, "....csv")
 
#paste results from one predictor under those of models made with two or more predictors
  result.Acartia <- resultaten.Acartia
  coeff.all.1   <- as.matrix(coeff.all.1)
  coeff.all.2   <- as.matrix(coeff.all.2)
  coeff.all     <- rbind(coeff.all.2,coeff.all.1)
  coeff.all     <- as.data.frame(coeff.all)
  coeff.all     <- as.matrix(coeff.all)
  colnames(coeff.all) <- c("Intercept","CHFLa","TEMP","SAL","PCB","Diatomratio")
  
  result.Acartia <- cbind(result.Acartia, coeff.all)  
  colnames(result.Acartia) <- c("Forms","AICs","AICcs","Intercept","CHFLa","TEMP","SAL","PCB","Diatomratio")

  
  #calculate weights per combination of predictors and store in a matrix
  result.Acartia <- result.Acartia[order(as.numeric(result.Acartia[,"AICcs"])),]         #sort according to aiccs
  tussenstappen <- NULL
  
  for (k in 1:nrow(result.Acartia))
  {
    delta.i <- as.numeric(result.Acartia[k,"AICcs"])- as.numeric(result.Acartia[1,"AICcs"])
    tussenstap <- exp(-delta.i/2)
    tussenstappen <-c(tussenstappen,tussenstap)
  }
  
  som.tussenstappen<-sum(tussenstappen)
  weigth.i <- tussenstappen/som.tussenstappen
  weigth.i <- weigth.i
  cum.weigth.i<-cumsum(weigth.i)
  
  #select those combinations where cumulated weight remains under 0.95, and paste the original weight of the models considering all possible combinations.
  selected.models <- which(cum.weigth.i<=0.95)   
  result.Acartia.sel <- result.Acartia[c(selected.models),]
  result.Acartia.sel <- cbind(result.Acartia.sel,weigth.i[selected.models])
  colnames(result.Acartia.sel)<-c("forms","AICs","AICcs","Intercept","CHFLa.org","TEMP.org","SAL.org","PCB.org","Diatomratio.org","weights.org")
  
  tussenstappen.sel <- NULL
  tussenstap.sel <- NULL
  
  
  
  #recalculate the weights only taking into account the selected models and paste this weight
  k <- NULL
  
  for (k in 1:nrow(result.Acartia.sel))
  {
    delta.i.sel <- as.numeric(result.Acartia.sel[k,"AICcs"])- as.numeric(result.Acartia.sel[1,"AICcs"])
    tussenstap.sel <- exp(-delta.i.sel/2)
    tussenstappen.sel <-c(tussenstappen.sel,tussenstap.sel)
  }
  
  som.tussenstappen.sel<-sum(tussenstappen.sel)
  weigth.i.sel <- tussenstappen.sel/som.tussenstappen.sel
 
  
  weigth.i.sel <- weigth.i.sel
  result.Acartia.sel <- cbind(result.Acartia.sel,weigth.i.sel)
  colnames(result.Acartia.sel)<-c("forms","AICs","AICcs","Intercept","CHFLa.org","TEMP.org","SAL.org","PCB.org","Diatomratio.org","weights.org","weights.sel")
  
  result.Acartia.sel <- result.Acartia.sel[,-10]
  
  #add matrix with weights per predictor and store in a matrix
  matrix.weights <- matrix(data=0,nrow=nrow(result.Acartia.sel),ncol=6)
  
  i <- NULL
  j <- NULL
  
  for (i in 1:6)
  {
    for (j in 1:nrow(result.Acartia.sel))
    {
      if (as.numeric(result.Acartia.sel[j,i+3])!= 0)
      { matrix.weights[j,i]<- as.numeric(result.Acartia.sel[j,"weights.sel"])
      }else{
        matrix.weights[j,i] <- 0
      }
    }
  }
  
  result.Acartia.sel <- cbind(result.Acartia.sel,matrix.weights[,1], matrix.weights[,2],matrix.weights[,3],matrix.weights[,4],matrix.weights[,5],matrix.weights[,6])
  colnames(result.Acartia.sel)<-c("forms","AICs","AICcs","Intercept.coef","CHFLa.coef","TEMP.coef","SAL.coef","PCB.coef","Diatomratio.coef","weights","Intercept.w","CHFLa.w","TEMP.w","SAL.w","PCB.w","Diatomratio.w")
  
  #calculate weighted.average and store in a matrix
  Intercept.weighted.average   <- as.numeric(result.Acartia.sel[,"Intercept.coef"])* as.numeric(result.Acartia.sel[,"weights"])
  CHFLa.weighted.average       <- as.numeric(result.Acartia.sel[,"CHFLa.coef"])* as.numeric(result.Acartia.sel[,"weights"])
  TEMP.weighted.average        <- as.numeric(result.Acartia.sel[,"TEMP.coef"])* as.numeric(result.Acartia.sel[,"weights"])
  SAL.weighted.average         <- as.numeric(result.Acartia.sel[,"SAL.coef"])* as.numeric(result.Acartia.sel[,"weights"])
  PCB.weighted.average         <- as.numeric(result.Acartia.sel[,"PCB.coef"])* as.numeric(result.Acartia.sel[,"weights"])
  Diatomratio.weighted.average <- as.numeric(result.Acartia.sel[,"Diatomratio.coef"])* as.numeric(result.Acartia.sel[,"weights"])

  
  result.Acartia.sel     <- cbind(result.Acartia.sel,Intercept.weighted.average,CHFLa.weighted.average,TEMP.weighted.average, SAL.weighted.average, PCB.weighted.average, Diatomratio.weighted.average)
  
  #add line with sums of colums
  extra.line <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,sum(abs(as.numeric(result.Acartia.sel[,"Intercept.w"]))),sum(abs(as.numeric(result.Acartia.sel[,"CHFLa.w"]))),sum(abs(as.numeric(result.Acartia.sel[,"TEMP.w"]))),sum(abs(as.numeric(result.Acartia.sel[,"SAL.w"]))),sum(abs(as.numeric(result.Acartia.sel[,"PCB.w"]))),sum(abs(as.numeric(result.Acartia.sel[,"Diatomratio.w"]))),sum(abs(as.numeric(result.Acartia.sel[,"Intercept.weighted.average"]))),sum(abs(as.numeric(result.Acartia.sel[,"CHFLa.weighted.average"]))),sum(abs(as.numeric(result.Acartia.sel[,"TEMP.weighted.average"]))),sum(abs(as.numeric(result.Acartia.sel[,"SAL.weighted.average"]))),sum(abs(as.numeric(result.Acartia.sel[,"PCB.weighted.average"]))),sum(abs(as.numeric(result.Acartia.sel[,"Diatomratio.weighted.average"]))))
  result.Acartia.sel <- rbind(result.Acartia.sel,extra.line)
  
  #calculate model average
  Intercept.model.averaged <- as.numeric(result.Acartia.sel["extra.line","Intercept.weighted.average"])/as.numeric(result.Acartia.sel["extra.line","Intercept.w"])
  CHFLa.model.averaged     <- as.numeric(result.Acartia.sel["extra.line","CHFLa.weighted.average"])/as.numeric(result.Acartia.sel["extra.line","CHFLa.w"])
  TEMP.model.averaged      <- as.numeric(result.Acartia.sel["extra.line","TEMP.weighted.average"])/as.numeric(result.Acartia.sel["extra.line","TEMP.w"])
  SAL.model.averaged       <- as.numeric(result.Acartia.sel["extra.line","SAL.weighted.average"])/as.numeric(result.Acartia.sel["extra.line","SAL.w"])
  PCB.model.averaged       <- as.numeric(result.Acartia.sel["extra.line","PCB.weighted.average"])/as.numeric(result.Acartia.sel["extra.line","PCB.w"])
  Diatomratio.model.averaged   <- as.numeric(result.Acartia.sel["extra.line","Diatomratio.weighted.average"])/as.numeric(result.Acartia.sel["extra.line","Diatomratio.w"])
  

  model.average <- matrix(data=c(Intercept.model.averaged,CHFLa.model.averaged,TEMP.model.averaged, SAL.model.averaged, PCB.model.averaged, Diatomratio.model.averaged),nrow=1,ncol=6)
  colnames(model.average) <-  c("Intercept.model.averaged","CHFLa.model.averaged","TEMP.model.averaged","SAL.model.averaged", "PCB.model.averaged","Diatomratio.model.averaged")
  
  #write a csv file of output
  result.Acartia.sel <- result.Acartia.sel[,-1]
  result.Acartia.sel <- as.matrix(result.Acartia.sel)
  
  formules <- result.Acartia[c(selected.models),"Forms"]
  
  #write.csv(formules,"....csv")
  #write.csv(result.Acartia.sel,"....csv") 
  #write.csv(model.average,"....csv") 