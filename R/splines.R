#' Mean curve estimation by fitting penalized splines to the data
#' 
#' @param elemID a vector of element ID
#' @param tVec a vector of observed time domain points
#' @param yVec a vector of response points
#' @param splineObj a spline object constructed on the time domain
#' @param lambda tuning parameter
#' @return a list of coefficients, residuals, and mean function
fitMeanCurve = function(obsData, splineObj, lambda){
    multiElement = "elemID" %in% colnames(obsData)
    measureError = "sigma" %in% colnames(obsData)
    
    tVec = obsData[,"obsT"]
    yVec = obsData[,"obsY"]
    if(measureError){
        ySigma = obsData[,"sigma"]
    }else{
        ySigma = NULL
    }
    
    if(!multiElement){
        meanModel = fitMeanCurveSingle(tVec, yVec, splineObj, 
                                       lambda, ySigma)
        meanModel = list(modelList = meanModel,
                         elemLevels = NULL)
    }else{
        elemID = obsData[,"elemID"]
        meanModel = fitMeanCurveMultiple(tVec, yVec, splineObj, 
                                         lambda, elemID, ySigma)
    }
    return(meanModel)
}


subtractMeanCurve = function(meanModel, inputData){
    outputData = inputData
    if(is.null(meanModel$elemLevels)){
        obsY = outputData[, "obsY"]
        obsT = outputData[, "obsT"]
        obsFitted = meanModel$modelList$fitSeq(obsT)
        outputData[selR, "obsY"] = obsY - obsFitted
    }else{
        eI = 0
        elemID = outputData[,"elemID"]
        for(bI in meanModel$elemLevels){
            eI = eI + 1
            selR = elemID == bI
            obsY = outputData[selR, "obsY"]
            obsT = outputData[selR, "obsT"]
            obsFitted = meanModel$modelList[[eI]]$meanFunction(obsT)
            outputData[selR, "obsY"] = obsY - obsFitted
        }
    }
    return(outputData)
}


plotMeanModel = function(meanModel, addPlot = FALSE){
    elemLevels = meanModel$elemLevels
    numPoints = 200

    plotData = data.frame()
    eI = 0
    for(eName in elemLevels){
        eI = eI + 1
        tmin = meanModel$modelList[[eI]]$tmin
        tmax = meanModel$modelList[[eI]]$tmax
        tSeq = seq(tmin, tmax, length.out = numPoints)
        ySeq = meanModel$modelList[[eI]]$meanFunction(tSeq)
        tmp = data.frame(eName, tSeq, ySeq)
        plotData = rbind(plotData, tmp)
    }
    colnames(plotData) = c("elemID", "obsT", "obsY")
    plotData$elemID = factor(plotData$elemID, levels = meanModel$elemLevels)
    
    if(addPlot){
        q = geom_line(aes(obsT, obsY), 
                      data = plotData, color = "black") 
    }else{
        q = ggplot(plotData, aes(obsT, obsY, color = elemID)) +
            geom_line()
    }
    return(q)
}



fitMeanCurveMultiple = function(tVec, yVec, splineObj, lambda,
                                elemID = NULL, ySigma = NULL){
    elemLevels = levels(elemID)
    modelList = list()
    numElem = length(elemLevels)
    if(length(lambda) == 1){
        lambda = rep(lambda, numElem)
    }else if(length(lambda) != numElem){
        stop("Incorrect dimension of lambda")
    }
    k = 0
    for(e in elemLevels){
        k = k + 1
        selR = (elemID == e)
        tVecSub = tVec[selR]
        yVecSub = yVec[selR]
        ySigmaSub = NULL
        if(!is.null(ySigma))
            ySigmaSub = ySigma[selR]
        modelS = fitMeanCurveSingle(tVecSub, yVecSub,  
                                    splineObj, lambda[k], ySigmaSub)
        modelList = c(modelList, list(modelS)) 
    }
    meanModel = list(elemLevels = elemLevels,
                     modelList = modelList)
    return(meanModel)    
}


fitMeanCurveSingle = function(tVec, yVec, 
                              splineObj, lambda,
                              ySigma = NULL){
    tmin = splineObj$getTMin()
    tmax = splineObj$getTMax()
    Omega = splineObj$get_Omega()
    Bmat = splineObj$evalSpline(tVec)
    if(!is.null(ySigma)){
        # Weight row by measurement error
        BmatW = diag(ySigma) %*% Bmat
        yVecW = ySigma * yVec
    }else{
        yVecW = yVec
        BmatW = Bmat
    }
    
    #compute coef
    tmp = t(BmatW) %*% BmatW + lambda * Omega
    beta = solve(tmp, t(BmatW) %*% yVecW)
    
    #compute fitted values
    residualVec = yVec - Bmat %*% beta
    residualVec = as.vector(residualVec)
    tSeq = seq(tmin, tmax, length.out = 200)
    ySeq = splineObj$evalSpline(tSeq) %*% beta
    fitSeq = approxfun(tSeq, ySeq)
    model = list(tmin = tmin, 
                 tmax = tmax,
                 beta = beta, 
                 meanFunction = fitSeq,
                 residualVec = residualVec)
    return(model)
}



