function write_final_header(fileName::String, rinfo::ReReturnInfo, prefixNames::Vector{String}, sep::String)
  open(fileName, "w") do f
    for name in prefixNames
      write(f, "$name$sep")
	  end
    write(f, "restart$(sep)")
    write(f, "popnSize$(sep)")
    write(f, "found$(sep)")
    write(f, "evals$(sep)")
    write(f, "bestRestart$(sep)")
    write(f, "bestFit$(sep)")
    write(f, "bestChr\n")
  end
end

function write_final(rinfo::ReReturnInfo; prefixValues = [], prefixNames = [], initialize = false,
	                 fileName = "results_withRestart", path = "", sep = ",", ext = ".txt")
  ext = (sep == ",") ? ".csv" : ext
  fileName = "$(path)/$fileName$ext"

  if initialize
    write_final_header(fileName, rinfo, prefixNames, sep)
  end

  state = last(rinfo.state)
  popnSize = last(rinfo.popnSize)
  totalEvals = last(rinfo.totalEvals)
  restart = length(rinfo.state) - 1
  foundSol = found(state)
  (bestChr, bestFit, bestRestart) = rinfo.best

  open(fileName, "a") do f
  	for value in prefixValues
  	  write(f, "$value$sep")
  	end
    write(f, "$restart$(sep)")
    write(f, "$popnSize$(sep)")
    write(f, "$foundSol$(sep)")
    write(f, "$totalEvals$(sep)")
    write(f, "$bestRestart$(sep)")
    write(f, "$bestFit$(sep)")
    write(f, "$bestChr\n")
  end
end

function write_final_header(fileName::String, rinfo::ReturnInfo, prefixNames::Vector{String}, sep::String)
  open(fileName, "w") do f
    for name in prefixNames
      write(f, "$name$sep")
	  end
    write(f, "popnSize$(sep)")
    write(f, "found$(sep)")
    write(f, "evals$(sep)")
    write(f, "bestFit$(sep)")
    write(f, "bestChr\n")
  end
end

function write_final(rinfo::ReturnInfo; prefixValues = [], prefixNames = [], initialize = false,
	                 fileName = "results_withRestart", path = "", sep = "\t", ext = ".txt")
  ext = (sep == ",") ? ".csv" : ext
  fileName = "$(path)/$fileName$ext"

  if initialize
    write_final_header(fileName, rinfo, prefixNames, sep)
  end

  state = rinfo.state
  popnSize = popnsize(state)
  totalEvals = evals(state)
  foundSol = found(state)
  (chr, fit) = best(state)

  open(fileName, "a") do f
  	for value in prefixValues
  	  write(f, "$value$sep")
  	end
    write(f, "$popnSize$(sep)")
    write(f, "$foundSol$(sep)")
    write(f, "$totalEvals$(sep)")
    write(f, "$fit$(sep)")
    write(f, "$chr\n")
  end
end
