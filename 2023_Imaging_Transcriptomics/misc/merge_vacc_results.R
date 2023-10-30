  args = commandArgs(trailingOnly = TRUE)
  sim_path=args[1]
  sim_n=as.numeric(args[2])

  remove.orig=T
  file_names<-sapply(1:sim_n, function(i) sprintf('%s/sim_%s_res.csv',sim_path,i))
  print(file_names)
  # check if all files exist
  test_exist=file.exists(file_names)
  if (all(test_exist)){
    theList = lapply(file_names, read.csv, stringsAsFactors=F)
    theResult = do.call(rbind,theList)

  
  write.csv(theResult, sprintf('%s/Res_%s_sim%s.csv',sim_path,basename(sim_path),length(theList)), row.names = F)
  if (remove.orig){
    tmp=lapply(file_names, file.remove)
  }
  } else {
    # find missing file
    missing_file=file_names[!test_exist]
    
    # get the file name
    missing_file=basename(missing_file)

    # extract number of missing_file from file name
    missing_num=as.numeric(gsub('.*sim_(.*)_res.csv','\\1',missing_file))
    # list missing files
    cat('missing files: \n')
    cat(missing_file)
  }
