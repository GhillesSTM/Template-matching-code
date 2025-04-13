function  Maxima =find_maxima(datafile, max, size, threshold)

      hLocalMax = vision.LocalMaximaFinder;% find local maxima
      hLocalMax.MaximumNumLocalMaxima = max;%Maximum number of maxima to find
      hLocalMax.NeighborhoodSize = [size size];%Specify the size of the neighborhood around the maxima
      hLocalMax.Threshold = threshold;%Value that all maxima should match or exceed

      Maxima = step(hLocalMax,datafile)%Find local maxima in input image