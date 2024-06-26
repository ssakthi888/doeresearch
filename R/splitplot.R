#' @title Analysis of Split plot design
#' @description The function gives ANOVA, R-square of the model, normality testing of residuals, SEd (standard error of difference), interpretation of ANOVA results and multiple comparison test for means.
#' @param data dependent variables
#' @param block vector containing replications
#' @param main.plot vector containing main-plot levels
#' @param sub.plot vector containing sub-plot levels
#' @param mean.comparison.test 0 for no test, 1 for LSD test, 2 for Dunccan test and 3 for HSD test
#' @param round.digits enter number of decimal digits to be displayed (0,1,2,3)
#' @param alpha select 0 for mean value only or 1 for mean value with mean.comparison.test
#' @importFrom agricolae LSD.test HSD.test duncan.test
#' @importFrom stats anova lm shapiro.test pf
#' @importFrom dplyr
#' @importFrom tidyr
#' @return ANOVA, interpretation of ANOVA, R-square, normality test result, SEd and multiple comparison test result
#' @export

splitplot <- function(data, block, main.plot, sub.plot, mean.comparison.test,round.digits,alpha) {
  if (alpha == 0) {
  if (mean.comparison.test == 0) {
    sak888 <- split2(data, block, main.plot, sub.plot, mean.comparison.test,round.digits)
  } else {
    results <- split2(data, block, main.plot, sub.plot, mean.comparison.test,round.digits)
    result2 <- data.frame()
    
    sak777 <- list()
    
    for (name in names(results)) {
      sak77 <- results[[name]][[1]][[14]]
      sak777 <- c(sak777, list(sak77))
    }
    
    sak1234 <- list()
    sak4321 <- list()
    for (name in names(results)) {
      if (mean.comparison.test == 1) {
        m_testAB <- results[[name]][[1]][[13]][[2]]
        main_test<- results[[name]][[1]][[9]][[2]]
        sub_test<- results[[name]][[1]][[11]][[2]]
        } else {
          m_testAB <- results[[name]][[1]][[13]][[3]]
          main_test<- results[[name]][[1]][[9]][[3]]
          sub_test<- results[[name]][[1]][[11]][[3]]
          
          }
      
      # Get existing row names as a new column
      existing_row_names <- rownames(m_testAB)
      m_testAB <- cbind(TRT = existing_row_names, m_testAB)
      
      # Set numeric row names
      rownames(m_testAB) <- 1:nrow(m_testAB)
      
      # Separate TRT into M and S columns based on ":"
      m_testAB <- m_testAB %>% separate(TRT, into = c("M", "S"), sep = ":")
      
      # Rename dependent.var to value and round it to 2 decimal places
      m_testAB <- m_testAB %>% rename(value = dependent.var)
      
      m_testAB$value <- round(m_testAB$value, round.digits)
      
      # Arrange by S and then M
      m_testAB <- m_testAB %>% arrange(S) %>% arrange(M)
      
      # Select relevant columns
      m_testAB <- m_testAB[, 1:3]
      
      # Pivot wider based on S, rename Treatments column
      result0 <- m_testAB %>%
        pivot_wider(names_from = S, values_from = value) %>%
        rename(Treatments = M)
      
      # Step 1: Compute mean for each numeric column and add a "Mean" row
      mean_row <- result0 %>%
        summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
        mutate(Treatments = "Mean") %>%
        mutate(across(where(is.numeric), ~ round(., round.digits)))
      
      # Step 2: Compute mean for each row across numeric columns and add a "Mean" column
      mean_col <- result0 %>%
        rowwise() %>%
        mutate(Mean = mean(c_across(where(is.numeric)), na.rm = TRUE)) %>%
        ungroup() %>%  # ungroup after rowwise operation
        mutate(across(where(is.numeric), ~ round(., round.digits)))
      
      # Combine mean_row and mean_col into the final result
      result <- bind_rows(mean_col, mean_row)
      
      # Step 3: Calculate the total mean (mean of mean row/mean column)
      total_mean_value <- mean(c(
        mean_row %>% select(where(is.numeric)) %>% unlist(), 
        mean_col %>% select(Mean) %>% unlist()
      ), na.rm = TRUE)
      
      # Step 4: Replace NA in the "Mean" column of the "Mean" row with total_mean_value
      result <- result %>%
        mutate(Mean = ifelse(Treatments == "Mean", round(total_mean_value, round.digits), Mean))

      # Convert all columns in result to characters
      result[] <- lapply(result, as.character)
       result5 <- result[]
       sak1234 = c(sak1234,list(result5))
       
       
       #main
       write.csv(main_test,"temp999.csv"); main_test=read.csv("temp999.csv"); file.remove("temp999.csv")
       main_test=main_test %>% arrange(X)  %>% rename(Main=dependent.var) %>% rename(Treatments=X)
       main_test$Main <- round(main_test$Main, round.digits)
       main_test=main_test %>% unite("Mean",Main:groups,sep="",na.rm=TRUE,remove = TRUE)
       
       
       #sub
       write.csv(sub_test,"temp999.csv"); sub_test=read.csv("temp999.csv"); file.remove("temp999.csv")
       sub_test=sub_test %>% arrange(X)  %>% rename(Sub=dependent.var) %>% rename(Treatments=X)
       sub_test$Sub <- round(sub_test$Sub, round.digits)
       sub_test=sub_test %>% unite("Mean",Sub:groups,sep="",na.rm=TRUE,remove = TRUE)
       sub_test=as.data.frame(t(sub_test),stringsAsFactors = FALSE)
       col_names <- make.unique(as.character(sub_test[1,]))
       colnames(sub_test) <- col_names
       sub_test <- sub_test[-1, ]
       write.csv(sub_test,"temp999.csv"); sub_test=read.csv("temp999.csv"); file.remove("temp999.csv")
       sub_test=sub_test %>% rename(Treatments=X)
       sub_test$Mean <- NA
       
       result_test1=merge(result0,main_test, by = "Treatments", all = TRUE)
       result_test1=rbind(result_test1,sub_test)
       result_test1 <- result_test1 %>%
         mutate(Mean = ifelse(Treatments == "Mean", round(total_mean_value, round.digits), Mean))
       
       # Convert all columns in result to characters
       result_test1[] <- lapply(result_test1, as.character)
       result_test <- result_test1[]
       sak4321 = c(sak4321,list(result_test))
    }

    sak888 <- list(results, sak1234,sak777) 
    
  }
  } else 
    if (alpha == 1){  if (mean.comparison.test == 0) {
    sak888 <- split2(data, block, main.plot, sub.plot, mean.comparison.test,round.digits)
  } else {
    results <- split2(data, block, main.plot, sub.plot, mean.comparison.test,round.digits)
    result2 <- data.frame()
    
    sak777 <- list()
    
    for (name in names(results)) {
      sak77 <- results[[name]][[1]][[14]]
      sak777 <- c(sak777, list(sak77))
    }
    
    sak1234 <- list()
    sak4321 <- list()
    for (name in names(results)) {
      if (mean.comparison.test == 1) {
        m_testAB <- results[[name]][[1]][[13]][[2]]
        main_test<- results[[name]][[1]][[9]][[2]]
        sub_test<- results[[name]][[1]][[11]][[2]]
      } else {
        m_testAB <- results[[name]][[1]][[13]][[3]]
        main_test<- results[[name]][[1]][[9]][[3]]
        sub_test<- results[[name]][[1]][[11]][[3]]
        
      }
      
      # Get existing row names as a new column
      existing_row_names <- rownames(m_testAB)
      m_testAB <- cbind(TRT = existing_row_names, m_testAB)
      
      # Set numeric row names
      rownames(m_testAB) <- 1:nrow(m_testAB)
      
      # Separate TRT into M and S columns based on ":"
      m_testAB <- m_testAB %>% separate(TRT, into = c("M", "S"), sep = ":")
      
      # Rename dependent.var to value and round it to 2 decimal places
      m_testAB <- m_testAB %>% rename(value = dependent.var)
      
      m_testAB$value <- round(m_testAB$value, round.digits)
      
      # Arrange by S and then M
      m_testAB <- m_testAB %>% arrange(S) %>% arrange(M)
      
      # Select relevant columns
      m_testAB <- m_testAB[, 1:3]
      
      # Pivot wider based on S, rename Treatments column
      result0 <- m_testAB %>%
        pivot_wider(names_from = S, values_from = value) %>%
        rename(Treatments = M)
      
      # Step 1: Compute mean for each numeric column and add a "Mean" row
      mean_row <- result0 %>%
        summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
        mutate(Treatments = "Mean") %>%
        mutate(across(where(is.numeric), ~ round(., round.digits)))
      
      # Step 2: Compute mean for each row across numeric columns and add a "Mean" column
      mean_col <- result0 %>%
        rowwise() %>%
        mutate(Mean = mean(c_across(where(is.numeric)), na.rm = TRUE)) %>%
        ungroup() %>%  # ungroup after rowwise operation
        mutate(across(where(is.numeric), ~ round(., round.digits)))
      
      # Combine mean_row and mean_col into the final result
      result <- bind_rows(mean_col, mean_row)
      
      # Step 3: Calculate the total mean (mean of mean row/mean column)
      total_mean_value <- mean(c(
        mean_row %>% select(where(is.numeric)) %>% unlist(), 
        mean_col %>% select(Mean) %>% unlist()
      ), na.rm = TRUE)
      
      # Step 4: Replace NA in the "Mean" column of the "Mean" row with total_mean_value
      result <- result %>%
        mutate(Mean = ifelse(Treatments == "Mean", round(total_mean_value, round.digits), Mean))
      
      # Convert all columns in result to characters
      result[] <- lapply(result, as.character)
      result5 <- result[]
      sak1234 = c(sak1234,list(result5))
      
      
      #main
      write.csv(main_test,"temp999.csv"); main_test=read.csv("temp999.csv"); file.remove("temp999.csv")
      main_test=main_test %>% arrange(X)  %>% rename(Main=dependent.var) %>% rename(Treatments=X)
      main_test$Main <- round(main_test$Main, round.digits)
      main_test=main_test %>% unite("Mean",Main:groups,sep="",na.rm=TRUE,remove = TRUE)
      
      
      #sub
      write.csv(sub_test,"temp999.csv"); sub_test=read.csv("temp999.csv"); file.remove("temp999.csv")
      sub_test=sub_test %>% arrange(X)  %>% rename(Sub=dependent.var) %>% rename(Treatments=X)
      sub_test$Sub <- round(sub_test$Sub, round.digits)
      sub_test=sub_test %>% unite("Mean",Sub:groups,sep="",na.rm=TRUE,remove = TRUE)
      sub_test=as.data.frame(t(sub_test),stringsAsFactors = FALSE)
      col_names <- make.unique(as.character(sub_test[1,]))
      colnames(sub_test) <- col_names
      sub_test <- sub_test[-1, ]
      write.csv(sub_test,"temp999.csv"); sub_test=read.csv("temp999.csv"); file.remove("temp999.csv")
      sub_test=sub_test %>% rename(Treatments=X)
      sub_test$Mean <- NA
      
      result_test1=merge(result0,main_test, by = "Treatments", all = TRUE)
      result_test1=rbind(result_test1,sub_test)
      result_test1 <- result_test1 %>%
        mutate(Mean = ifelse(Treatments == "Mean", round(total_mean_value, round.digits), Mean))
      
      # Convert all columns in result to characters
      result_test1[] <- lapply(result_test1, as.character)
      result_test <- result_test1[]
      sak4321 = c(sak4321,list(result_test))
    }
    sak888 <- list(results, sak4321,sak777) 
    
  }} else {
    stop("alpha must be 0 or 1")
  }
  return(sak888)
}

split2<-function(data,block,main.plot,sub.plot,mean.comparison.test,round.digits){
  split1<-function(dependent.var,block,main.plot,sub.plot,mean.comparison.test,round.digits,alpha){
    dependent.var<-as.numeric(dependent.var)
    block<-as.factor(block)
    main.plot<-as.factor(main.plot)
    sub.plot<-as.factor(sub.plot)
    model<-lm(dependent.var~block+block:main.plot+main.plot*sub.plot)
    n.test<-shapiro.test(model$residuals)
    n.result<-n.test$p.value
    if (n.result > 0.05){
      n.result<-"Normality assumption is not violated"
    } else {
      n.result<-"Normality assumption is violated"
    }
    summary.model<-summary(model)
    R.square<-paste("R Square",round(summary.model$r.squared,digits = 3))
    anova1<-anova(model)
    a<-anova1[3,]
    anova1[3,]<-anova1[4,]
    anova1[4,]<-a
    b<-rownames(anova1)[3]
    rownames(anova1)[3]<-"Ea"
    rownames(anova1)[4]<-b
    row.names(anova1)[6]<-"Eb"
    anova1[1,4]<-anova1[1,3]/anova1[3,3]
    anova1[2,4]<-anova1[2,3]/anova1[3,3]
    anova1[3,4]<-NA
    anova1[3,5]<-NA
    anova1[1,5]<-pf(as.numeric(anova1[1,4]),as.numeric(anova1[1,1]),as.numeric(anova1[3,1]), lower.tail = FALSE)
    anova1[2,5]<-pf(as.numeric(anova1[2,4]),as.numeric(anova1[2,1]),as.numeric(anova1[3,1]), lower.tail = FALSE)
    anova1[4,5]<-pf(as.numeric(anova1[4,4]),as.numeric(anova1[4,1]),as.numeric(anova1[6,1]), lower.tail = FALSE)
    anova1[5,5]<-pf(as.numeric(anova1[5,4]),as.numeric(anova1[5,1]),as.numeric(anova1[6,1]), lower.tail = FALSE)
    
    # Coefficient of Variation
    CVa <- sqrt(as.numeric(anova1[3, 3])) / mean(dependent.var) * 100
    CVa <- round(CVa,3)
    CVb <- sqrt(as.numeric(anova1[6, 3])) / mean(dependent.var) * 100
    CVb <- round(CVb,3)
    CV<-paste("CV(main):",CVa,", CV(sub):",CVb, ", CV(M x S):",CVb,", CV(S x M):",CVb)
    
    
    if (anova1[2,5]>0.05){
      multiple.comparison.A<-"All the main plot factor level means are same so dont go for any multiple comparison test"
    } else {
      multiple.comparison.A<-"The means of one or more levels of main plot factor are not same, so go for multiple comparison test"
    }
    if (anova1[4,5]>0.05){
      multiple.comparison.B<-"All the sub plot factor factor level means are same so dont go for any multiple comparison test"
    } else {
      multiple.comparison.B<-"The means of one or more levels of sub plot factor are not same, so go for multiple comparison test"
    }
    if (anova1[5,5]>0.05){
      multiple.comparison.AB<-"All the interaction level means are same so dont go for any multiple comparison test"
    } else {
      multiple.comparison.AB<-"The means of one or more levels of interaction are not same, so go for multiple comparison test"
    }
    
    # Calculate SEd
    replications <- length(unique(block))
    n_main_plot <- length(unique(main.plot))
    n_sub_plot <- length(unique(sub.plot))
    
    # SEd for main plot
    MSa <- anova1["Ea", "Mean Sq"]
    SEd_main <- sqrt(2 * MSa / (replications * n_sub_plot))
    SEd_main <-round(SEd_main,3)
      
    # SEd for sub plot
    MSb <- anova1["Eb", "Mean Sq"]
    SEd_sub <- sqrt(2 * MSb / (replications * n_main_plot))
    SEd_sub <-round(SEd_sub,3)
    
    # SEd for interaction plot
    SEd_MxS <- sqrt(2 * MSb / replications)
    SEd_MxS <-round(SEd_MxS,3)
    
    # SEd for M at the same level of S
    SEd_SxM <- sqrt(2 * (MSa + (n_sub_plot - 1) * MSb) / (replications * n_sub_plot))
    SEd_SxM <-round(SEd_SxM,3)
    
    # Combine SEd values into a data frame
    SEd<-paste("SEd(main):",SEd_main,", SEd(sub):",SEd_sub, ", SEd(M x S):",SEd_MxS,", SEd(S x M):",SEd_SxM)
    
    # Degrees of freedom
    df_main <- anova1["Ea", "Df"]
    df_sub <- anova1["Eb", "Df"]
    df_MxS <- anova1["Eb", "Df"]
    
    # Critical t-values for 95% confidence (two-tailed)
    t_main <- qt(0.975, df_main)
    t_sub <- qt(0.975, df_sub)
    t_MxS <- qt(0.975, df_MxS)
    t_SxM <- ((t_main*MSa)+(t_sub*(n_sub_plot-1)*MSb))/(MSa+((n_sub_plot-1)*MSb))
    
    # Construct t_test data frame
    t_test = c(t_main, t_sub, t_MxS, t_SxM)
    
    # Construct CD 
    CD_main = round(t_main*SEd_main,3)
    CD_sub = round(t_sub*SEd_sub,3)
    CD_MxS = round(t_MxS*SEd_MxS,3)
    CD_SxM = round(t_SxM*SEd_SxM,3)
    
    # Combine CD values into a data frame
    CD<-paste("CD(main):",CD_main,", CD(sub):",CD_sub, ", CD(M x S):",CD_MxS,", CD(S x M):",CD_SxM)
   
    #dataframe
    CV1 = c(CVa, CVb, CVb, CVb)
    SEd1 = c(SEd_main, SEd_sub, SEd_MxS, SEd_SxM)
    t_test1 = c(t_main, t_sub, t_MxS, t_SxM)
    CD1 = c(CD_main, CD_sub, CD_MxS, CD_SxM)
    
    sak123 <- data.frame(Category = c("Main", "Sub", "M x S", "S x M"),
                         SEd=round(SEd1,round.digits),CD=round(CD1,round.digits))
    # Add significance indication
    significant <- ifelse(anova1[c(2, 4, 5, 5), 5] <= 0.05, "S", "NS")
    sak123$CD <- ifelse(significant == "NS", "NS", sak123$CD)
    sak123 <- sak123 %>% rename("CD (P=0.05)"=CD)
    
    sak123 <- as.data.frame(t(sak123))
    colnames(sak123) <- sak123[1, ]
    sak123<- sak123[-1, ]
    
    # Get existing row names as a new column
    existing_row_names1 <- rownames(sak123)
    sak123 <- cbind(LSD=existing_row_names1, sak123)
    # Create an empty column as the first column
   # "NA" <- rep("", nrow(sak123))
   # names(empty_col) <- ""  # Set an empty name to avoid displaying any column name
    # Add the empty column to ft2 data frame
    #sak123 <- cbind("NA", sak123)

    # Set numeric row names
    rownames(sak123) <- 1:nrow(sak123) 
    
    
    if (mean.comparison.test == 0){
      m.testA<-"No multiple comparison test selected"
      m.testB<-m.testA
      m.testAB<-m.testA
      sak123<-m.testA
    }
    if (mean.comparison.test == 1){
      m.test1A<-LSD.test(dependent.var,main.plot,anova1[3,1],anova1[3,3])
      m.testA<-list(m.test1A$statistics,m.test1A$groups)
      m.test1B<-LSD.test(dependent.var,sub.plot,anova1[6,1],anova1[6,3])
      m.testB<-list(m.test1B$statistics,m.test1B$groups)
      m.test1AB<-LSD.test(dependent.var,main.plot:sub.plot,anova1[6,1],anova1[6,3])
      m.testAB<-list(m.test1AB$statistics,m.test1AB$groups)
    }
    if (mean.comparison.test == 2){
      m.test1A<-duncan.test(dependent.var,main.plot,anova1[3,1],anova1[3,3])
      m.testA<-list(m.test1A$statistics,m.test1A$duncan,m.test1A$groups)
      m.test1B<-duncan.test(dependent.var,sub.plot,anova1[6,1],anova1[6,3])
      m.testB<-list(m.test1B$statistics,m.test1B$duncan,m.test1B$groups)
      m.test1AB<-duncan.test(dependent.var,main.plot:sub.plot,anova1[6,1],anova1[6,3])
      m.testAB<-list(m.test1AB$statistics,m.test1AB$duncan,m.test1AB$groups)
    }
    if (mean.comparison.test == 3){
      m.test1A<-HSD.test(dependent.var,main.plot,anova1[3,1],anova1[3,3])
      m.testA<-list(m.test1A$statistics,m.test1A$parameters,m.test1A$groups)
      m.test1B<-HSD.test(dependent.var,sub.plot,anova1[6,1],anova1[6,3])
      m.testB<-list(m.test1B$statistics,m.test1B$parameters,m.test1B$groups)
      m.test1AB<-HSD.test(dependent.var,main.plot:sub.plot,anova1[6,1],anova1[6,3])
      m.testAB<-list(m.test1AB$statistics,m.test1AB$parameters,m.test1AB$groups)
    }
    
    
    my.list<-list(anova1,CD,SEd,CV,R.square,n.test,n.result,multiple.comparison.A,m.testA,multiple.comparison.B,m.testB,multiple.comparison.AB,m.testAB,sak123)
    done<-list(my.list)
    return(done)
  }
  fiftn<-convert(data)
  colnumber<-ncol(data)
  output<-list()
  for (j in 1:colnumber){
    output[[j]]<-split1(fiftn[[j]],block,main.plot,sub.plot,mean.comparison.test,round.digits)
  }
  names(output)<-names(data)
  return(output)
}


convert<-function(data1){
  data1<- as.data.frame(sapply(data1, as.numeric))
  data1<-as.list(data1)
  return(data1)
}
