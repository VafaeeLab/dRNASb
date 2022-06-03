#returns a list of de top table for hour 2 against hour 0,
#             hour4Vshour0, hour8Vshour0, hour16Vshour0, hour24Vshour0
de_analysis_hourwise <- function(data, pheno, method = "limma"){
  if(method == "limma"){

    data <- data.norm

    group <- as.factor(pheno[colnames(data), "groups"])
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    rownames(design) <- colnames(data)
    fit <- limma::lmFit(data, design = design)
    cont.matrix <- limma::makeContrasts(WT.02_h - WT.00_h,
                                          WT.04_h - WT.00_h,
                                          WT.08_h - WT.00_h,
                                          WT.16_h - WT.00_h,
                                          WT.24_h - WT.00_h,
                                        levels = design)
    fit.cont <- limma::eBayes(limma::contrasts.fit(fit, cont.matrix))
    summa.fit <- limma::decideTests(fit.cont)

    obtain_de_per_coef <- function(fit.cont, coef = 1, lfc = 1, adjP = 0.05){
      wtdt <- limma::topTable(fit.cont, n = Inf, coef = coef)
      updt = with(wtdt, logFC > lfc & adj.P.Val < adjP)
      downdt = with(wtdt, logFC < lfc & adj.P.Val < adjP)
      return (wtdt)
    }

    DE <- list()
    for (n in 1:5){
      DE[[n]] <- obtain_de_per_coef(fit.cont, coef = n)
    }

    return (DE)

  }

}
