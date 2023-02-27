# Wrappers to call external commands using system2
# Currently implemented: cdo, nco

cdo <-
  function(command,
             args = "",
             input = "",
             options = "",
       tmpdir = "",
             output = "",
             stdout = "",
             noout = F) {
    if (args != "") {
      args <- paste0(",", args)
    }
    if (stdout != "") {
      stdout <- paste0(" > '", stdout, "'")
      noout <- T
    }
    if (input[1] != "") {
      for (i in seq_along(input)) {
        input[i] <- paste0("'", input[i], "'")
      }
      input <- paste(input, collapse = " ")
    }
    output0 <- output
    if (output != "") {
      output <- paste0("'", output, "'")
    } else if (!noout) {
   if(tmpdir == ""){
    tmpdir=getwd()
      output <- tempfile(tmpdir=tmpdir)
   } # if MG
      output0 <- output
    }
    argstr <- paste0(
      options, " ", command, args, " ", input, " ", output,
      " ", stdout
    )
    ret <- system2("cdo", args = argstr)
    if (ret != 0) {
      stop(paste("Failed (", ret, "): cdo", argstr))
    }
    return(output0)
  }

#--- Coded by Min-gyu (20.12.15)
nco <-
  function(cmd,
             args = "",
             input = "",
             options = "",
       tmpdir = "",
             output = "",
             stdout = "",
             noout = F) {
    if (args != "") {
      args <- paste0(args)
    }
    if (stdout != "") {
      stdout <- paste0(" > '", stdout, "'")
      noout <- T
    }
    if (input[1] != "") {
      for (i in seq_along(input)) {
        input[i] <- paste0("'", input[i], "'")
      }
      input <- paste(input, collapse = " ")
    }
    input <- noquote(input)
    output0 <- output
    if (output != "") {
      output <- paste0("'", output, "'")
    } else if (!noout) {
   if(tmpdir == ""){
    tmpdir=getwd()
      output <- tempfile(tmpdir=tmpdir)
      output <- paste0("'",output,"'")
   } # if MG
      output0 <- output
    }
    argstr <- paste0(
      " ",options, " ", args, " ", input, " ", output)
    print(noquote(paste0(cmd, argstr)))
    return(output0)
  }
