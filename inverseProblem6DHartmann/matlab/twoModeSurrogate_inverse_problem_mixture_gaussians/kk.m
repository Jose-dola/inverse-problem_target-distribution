er = 100000;
while true
  inverse_problem_gaussians_twoModeSurrogate;
  fi = fopen('bestChainElemntErrorBefore.double.bin','r');
  er = fread(fi, [1 1] , 'double');
  fclose(fi);
  if modelerr(bestChain) < er
    chosenChain = chain;
    er                    
    er = modelerr(bestChain)   
  end
  !mv bestChainElemntError.double.bin bestChainElemntErrorBefore.double.bin                           
end

