
findPeakHeightsFromIndices = function(doubletPeakIndices, peakCenters, rawData){
  peakParamDataFrame = t(sapply(doubletPeakIndices, findAverageDoubletHeight, peakCenters=peakCenters, rawData = rawData))
  colnames(peakParamDataFrame) = c("massRaw","massRaw1","massRaw2","baseline","nodeDev", "baselineSlope","gap","antinodeDifferenceRaw")
  peakParamDataFrame
}

findIndicesBelowThreshold = function(x, t) which(x<t, arr.ind=T)
findPeakStartIndices = function(i) i[c(1, which(diff(i)> 1, arr.ind=T)+1)]
findPeakEndIndices = function(i) i[c(which(diff(i)> 1, arr.ind=T), length(i))]
minimumWithinRange = function(i1, i2, x) min(x[min(i1,i2):max(i1,i2)])
findPeakHeights = function(i1, i2, x) mapply(minimumWithinRange, i1, i2, MoreArgs=list(x=x))
which.minWithinRange = function(i1, i2, x) which.min(x[min(i1,i2):max(i1,i2)])+min(i1,i2)-1
findPeakCenters = function(i1, i2, x) mapply(which.minWithinRange, i1, i2, MoreArgs=list(x=x))


findDoubletPeakIndices = function(diffTimes, heights){	
  nextHeight = c(heights[-1], 1e5)
  timeToPrev = c(1e4, diffTimes)
  timeToNext = c(diffTimes, 1e4)
  timeAfterNext = c( diffTimes[-1],1e4, 1e4)
  which(timeToNext > tDoubletGap[1] & timeToNext < tDoubletGap[2] & abs(heights-nextHeight) < maxHeightDiff)  
}

findAverageDoubletHeight = function(i, rawData, peakCenters, plotTheFit=FALSE){
  dataIndicesInRegion = max((peakCenters[i]-baselineSearchWidthSidePoints),1) : min(peakCenters[i+1]+baselineSearchWidthSidePoints, length(rawData))
  dataInRegion = rawData[dataIndicesInRegion]
  peak1Center = baselineSearchWidthSidePoints  # relative to region
  peak2Center = baselineSearchWidthSidePoints+(peakCenters[i+1]-peakCenters[i])  # relative to region
  
  # truncate it so it has evenly sized chunks
  chunkSize = floor(length(dataInRegion)/numBaselineChunks)
  dataInRegion = dataInRegion[1:(numBaselineChunks*chunkSize)]
  dataIndicesInRegion = dataIndicesInRegion[1:(numBaselineChunks*chunkSize)]
  chunkSDs = apply(matrix(dataInRegion, ncol=numBaselineChunks), 2,sd)
  
  preBaselineChunkIndex = which.min(chunkSDs[1:floor(peak1Center/chunkSize)])
  postBaselineChunkIndex = which.min(chunkSDs[ceiling(peak2Center/chunkSize):length(chunkSDs)])+ceiling(peak2Center/chunkSize)-1
  
  tipSidePoints = round(tipSidePointFractionOfPeakWidth*(peak2Center-peak1Center))
  baselineSegmentIndices = c( (preBaselineChunkIndex-1)*chunkSize + 1:chunkSize, 
                              (postBaselineChunkIndex-1)*chunkSize + 1:chunkSize )
  baselineSegments = rep(NA, length(dataInRegion))
  baselineSegments[baselineSegmentIndices] = dataInRegion[baselineSegmentIndices]
  t = 1:length(dataInRegion)
  baselineModel = lm(baselineSegments~t)
  baseline = predict(baselineModel, newdata=data.frame(t))
  
  dataInRegion = dataInRegion - baseline
  antinodeRegion1 = dataInRegion[ peak1Center + (-tipSidePoints):tipSidePoints+1 ]
  antinodeRegion2 = dataInRegion[ peak2Center + (-tipSidePoints):tipSidePoints+1 ]
  nodeRegion=dataInRegion[ peak1Center:peak2Center ]
  t = 1:(2*tipSidePoints+1)
  mod1 = lm(antinodeRegion1 ~ stats::poly(t,4))
  mod2 = lm(antinodeRegion2 ~ stats::poly(t,4))
  tripletPeakHeight = (min(predict(mod1)) + min(predict(mod2)))/2
  
  if (plotTheFit){
    plot(freqTime[dataIndicesInRegion], dataInRegion,xlab='Time (s)', ylab='buoy. mass (~pg)', cex=.7)
    lines(freqTime[dataIndicesInRegion], freqData_bandPassFiltered[dataIndicesInRegion], col='red')
    abline(h=0,col='grey',lty=2)
    points(freqTime[dataIndicesInRegion], baselineSegments-baseline, col='blue')
    lines(freqTime[t + peak1Center-tipSidePoints + min(dataIndicesInRegion) - 1], predict(mod1), col='green', lwd=2)
    lines(freqTime[t + peak2Center-tipSidePoints + min(dataIndicesInRegion) - 1], predict(mod2), col='green', lwd=2)
  }
  
  output_struct <-  c(-tripletPeakHeight, # mass
                      -min(predict(mod1)),# massRaw1
                      -min(predict(mod2)),  # massRaw2
                      mean(baseline[c(peak1Center,peak2Center)]), # baseline
                      max(nodeRegion), # nodeDev
                      coef(baselineModel)[2], #baselineSlope
                      (peakCenters[i+1]-peakCenters[i]), #gap
                      min(predict(mod2))-min(predict(mod1))) #antinodeDifferenceRaw
  
  return(output_struct)
}	

getPeaks = function(freqData_bandPassFiltered, rawData, detectionThreshold, 
                    tDoubletGap,  maxHeightDiff, baselineSearchWidthSidePoints, numBaselineChunks,
                    tipSidePointFractionOfPeakWidth, visualizeNum=5) {
  
  freqTime = (1:length(freqData_bandPassFiltered))/dataRate # seconds
  
  indicesBelowThreshold = findIndicesBelowThreshold(freqData_bandPassFiltered, t=detectionThreshold)
  peakStartIndices = findPeakStartIndices(indicesBelowThreshold)
  peakEndIndices = findPeakEndIndices(indicesBelowThreshold)
  peakWidth_indices = peakEndIndices - peakStartIndices
  peakHeights = findPeakHeights(peakStartIndices, peakEndIndices, freqData_bandPassFiltered)
  peakCenters = findPeakCenters(peakStartIndices, peakEndIndices, freqData_bandPassFiltered)	
  
  ## Sometimes when the cell passes the node, the signal won't return to baseline.
  ## We need to detect those cases (detected as a single peak) and split them into two peaks.
  # INPUTS: peakHeights, freqData_bandPassFiltered, peakEndIndices
  
  toInsert = list()
  j = 1
  for (i in 1:length(peakHeights)){
    #threshold at half peak height
    indicesBelowThreshold = which(freqData_bandPassFiltered[peakStartIndices[i]:peakEndIndices[i]] < 0.85*peakHeights[i])
    newPeakEndIndices = indicesBelowThreshold[which(diff(indicesBelowThreshold) > tDoubletGap[1])]
    
    #if no discontinuities, continue
    if(length(newPeakEndIndices) != 1) next
    
    newPeakStartIndices = indicesBelowThreshold[which(diff(indicesBelowThreshold) > tDoubletGap[1],arr.ind=T)+1]
    toInsert[[j]] = c(i, peakStartIndices[i]+newPeakStartIndices, peakStartIndices[i]+newPeakEndIndices)
    j=j+1
  }
  
  # insert one at at time from end to ensure we don't mess up indexing
  for (ins in rev(toInsert)){ 
    if (ins[1]==length(peakStartIndices)){ 
      peakStartIndices = c(peakStartIndices,ins[2])
    } else {
      peakStartIndices = c(peakStartIndices[1:ins[1]],ins[2], peakStartIndices[(ins[1]+1):length(peakStartIndices)])
    }
    if (ins[1]==1) {
      peakEndIndices =   c(ins[3], peakEndIndices)
    } else {
      peakEndIndices =   c(peakEndIndices[1:(ins[1]-1)], ins[3], peakEndIndices[ins[1]:length(peakEndIndices)])
    }
  }	
  
  # now re-estimate these	
  peakHeights = findPeakHeights(peakStartIndices, peakEndIndices, freqData_bandPassFiltered)
  peakCenters = findPeakCenters(peakStartIndices, peakEndIndices, freqData_bandPassFiltered)	
  
  rm(indicesBelowThreshold)
  
  ####
  # Everything below appears to just look at peakHeights and peakCenters? (and timeBetweenPeaks, doubletPeakIndices)
  
  
  # GATING AND REFINING ESTIMATION
  # find all peaks that meet these criteria:  
  #    2 peaks within tDoubletGap[2] datapoints of each other and >tDoubletGap[1] datapoints between them. 
  
  timeBetweenPeaks = diff(peakCenters)
  doubletPeakIndices = findDoubletPeakIndices(timeBetweenPeaks, peakHeights)  
  doubletPeakParams = findPeakHeightsFromIndices(doubletPeakIndices, peakCenters, rawData)
  
  # plot visualizeNum random peaks along with fits
  if (visualizeNum!=0){
    for (i in 1:visualizeNum){
      j = sample(doubletPeakIndices,1)
      findAverageDoubletHeight(j, rawData, peakCenters, plotTheFit=TRUE)
    }
  }
  
  if (length(doubletPeakIndices)==0) return(data.frame(time=NULL, mass=NULL)) # in case it's empty
  peakParamDataFrame =  data.frame(time = freqTime[peakCenters[doubletPeakIndices]]/3600,	doubletPeakParams)
  
  return (list(singletsDetected=length(peakHeights), peakParamDataFrame=peakParamDataFrame)) 
}