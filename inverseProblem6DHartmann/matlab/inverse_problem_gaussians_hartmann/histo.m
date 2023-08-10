classdef histo < handle
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %   histo structure (handle): 
  %     histo.first : first point of the histogram
  %     histo.last  : last point of the histogram
  %     histo.n : number of intervals in the histogram
  %     histo.L : total length (last point - first point)
  %     histo.l : interval length
  %     histo.x : vector of histogram points : [first point of the first interval, middle point 
  %               of the first interval, first point of the second interval (and last point of
  %               first interval), middle point of the second interval, first point of the
  %               third interval (and last point of the second interval), ...]
  %     histo.w : (before calling normalize) vector of intervals' weights : [number of observations 
  %               in the first interval, number of observations in the second interval, number of
  %               observations in the third interval, ......., number of observations that are out
  %               of the defined intervals (p<histo.first or p>histo.last)]
  %               OR
  %               (after calling normalize) vector of intervals' densities : [probability of the  
  %               first interval, probability of the second interval, ...., probability of the
  %               points out of the defined intervals]
  %     histo.lowest  : lowest observation
  %     histo.highest : highest observation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  properties
    first;
    last;
    n;
    x;
    L;
    l;
    w;
    lowest;
    highest;
  end
  methods
    function histo = histo(first,last,numberOfIntervals)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % param:
      %   first : first point of the histogram
      %   last  : last point of the histogram
      %   numberOfIntervals : number of intervals in the histogram
      % return:
      %   histo structure (handle) according to the parameters.
      %   OBSERVATION: variable 'w' of the histo structure is set as a 
      %   zero vector by this function.
      %   OBSEVATION: variable 'lowest' is set as the first point of the histogram
      %   OBSEVATION: variable 'highest' is set as the last point of the histogram
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      histo.first   = first;
      histo.lowest  = first;
      histo.last    = last;
      histo.highest = last;
      histo.n       = numberOfIntervals;
      histo.x       = linspace(first,last,numberOfIntervals*2+1);
      histo.L       = last-first;
      histo.l       = histo.L/numberOfIntervals;
      histo.w       = zeros(1,numberOfIntervals+1);
    end
    function pointToHisto(histo,p)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % include the p observation in the histogram represented 
      % by the histo structure (handle)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if p > histo.last
        index = histo.n+1;
        if p > histo.highest
          histo.highest = p;
        end
      elseif p < histo.first
        index = histo.n+1;
        if p < histo.lowest
          histo.lowest = p;
        end
      else
        index = floor((p-histo.first)/histo.l)+1;
      end
      histo.w(index) = histo.w(index)+1;
    end
    function resetObservations(histo)
      histo.w = zeros(1,histo.n+1);
      histo.lowest  = histo.first;
      histo.highest = histo.last;
    end
    function normalize(histo)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % calculate interval densities of the histo structure (handle)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      l   = histo.first - histo.lowest + histo.highest - histo.last;
      aux = histo.l*sum(histo.w(1:end-1)) + l*histo.w(end);
      histo.w = histo.w ./ aux;
    end
    function e = L2error(histo,trueDensity)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % calculate L2 norm between the density represented by
      % the histogram and the especified PDF. We assume that
      % the especified PDF has no density out of the intervals
      % that are especified by the histogram.
      % OBSERVATION: this method should be called after calling
      % the 'normalize' method.
      %
      % param:
      %   trueDensity : vector with the evaluations of the PDF 
      %   in the 'middle points' of the 'histo.x' variable 
      %   (indices n such that n%2 == 0). You can call the
      %   'intervalrepresentativePoints' method to get these points
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      l = histo.first - histo.lowest + histo.highest - histo.last;
      e = histo.l*(sum((trueDensity-histo.w(1:end-1)).^2)) + l*histo.w(end)*histo.w(end);
    end
    function x = intervalrepresentativePoints(histo)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % return the 'middle points' of the intervals
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      x = histo.x(mod(1:length(histo.x),2)==0);
    end
    function plot(histo,color)
      plot( histo.x(mod(1:length(histo.x),2)==0) , histo.w(1:end-1) , color );
    end
  end
end
