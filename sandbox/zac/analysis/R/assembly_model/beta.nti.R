beta.nti.calc.stegen = function (samp,reps,path.to.reps,beta.mntd.obs) {
  
  rand.beta.mntd.comp = array(c(-999),dim=c(nrow(samp),nrow(samp),length(reps))); print(dim(rand.beta.mntd.comp)); # array to hold randomizations
  
  for (i in 1:length(reps)) { ## compile randomizations into one matrix
    
    temp.rand = as.matrix(read.csv(paste(path.to.reps,reps[i],'.csv',sep=""),row.names=1,header=T)); ## read in each replicate of randomized beta.mntd and place it in the array
    colnames(temp.rand) = gsub("X","",colnames(temp.rand))
    #colnames(temp.rand) = gsub('X','',colnames(temp.rand));
    if (identical(colnames(temp.rand),rownames(temp.rand))==F) {
      
      print('Error: temp.rand col/row name mismatch')
      print(i)
      
      break()}; # end if statement
    
    rand.beta.mntd.comp[,,i] = temp.rand;
    print(c(i,colnames(beta.mntd.obs)[colnames(beta.mntd.obs) != colnames(temp.rand)])); # looking for any mismatches
    rm('temp.rand');
    
  }; 
  
  rm("i");
  
  print(rand.beta.mntd.comp[1:5,1:5,1])
  print(rand.beta.mntd.comp[1:5,1:5,2])
  print(dim(rand.beta.mntd.comp))
  
  beta.nti = matrix(c(NA),nrow=nrow(samp),ncol=nrow(samp));
  
  for (columns in 1:(nrow(samp)-1)) {
    for (rows in (columns+1):nrow(samp)) {
      
      rand.vals = rand.beta.mntd.comp[rows,columns,];
      beta.nti[rows,columns] = (beta.mntd.obs[rows,columns] - mean(rand.vals)) / sd(rand.vals);
      rm("rand.vals");
      
    };
  };
  
  rownames(beta.nti) = rownames(beta.mntd.obs); colnames(beta.nti) = colnames(beta.mntd.obs);
  return(beta.nti)
  
}; # end beta.nti.calc.stegen function
