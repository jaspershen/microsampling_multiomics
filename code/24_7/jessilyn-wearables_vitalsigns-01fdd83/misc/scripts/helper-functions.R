bootstrap.dataset = function(data, replace = FALSE){
  n = nrow(data)
  data[sample(n,replace = replace),]
}
correct.vars = function(data){
  vars = list("HSCRP","TBIL","IGM","TGL","UALB","AG","AST","LDL")
  vars = vars[vars %in% colnames(data)]
  for (var in vars){
    data[,var] = gsub("<", "", data[,var])
    data[,var] = as.numeric(data[,var])
  }
  data
}


