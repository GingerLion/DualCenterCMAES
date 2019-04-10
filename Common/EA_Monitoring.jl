# abstract type Monitor end
# abstract type RunInfo <: Monitor end
#
# struct NoMonitor <: Monitor end
#
# mutable struct EA_RunInfo <: RunInfo
#   data::Dict
# end
#
# const Storage = RunInfo           # depreciated
# const NoStorage = NoMonitor       # depreciated

RunInfo() = EA_RunInfo()
RunInfo(ri::RunInfo) = RunInfo()    # RunInfo can be extended with constants that can be copied - no constants for a general EA
NoMonitor(ri::RunInfo) = NoMonitor()

function monitor!(ns::NoMonitor, state::State) end

# ds = DetailsStore()
monitorable(ri::RunInfo) = true
monitorable(nm::NoMonitor) = false

newmonitor(monitoring::Bool) =  (monitoring ? RunInfo() : NoMonitor())
newmonitor(ri::Monitor) =  typeof(ri)(ri)

function reverse!(ri::RunInfo)
  for key in keys(ri.data)
    ri.data[key] = reverse(ri.data[key])
  end
  ri
end

reverse!(nm::NoMonitor) = nm

length(ri::RunInfo) = length(ri[:gen])

function collect!(ri::RunInfo)
  for key in keys(ri.data)
    ri.data[key] = collect(ri.data[key])
  end
  ri
end

collect!(nm::NoMonitor) = nm


function setupmonitor(detailed::Bool)
  detailed ? RunInfo() : NoMonitor()
end

function key_exists(dict::Dict, keyName)
  typeof(getkey(dict, keyName, false)) != Bool
end

function add_values!(ri::RunInfo; keyName = "w", values = nil())
  # ke = key_exists(ds.data, keyName)
  # println("key_exists = $(ke)")
  if key_exists(ri.data, keyName)
      vals = cons(values, ri.data[keyName])
  else
      vals = list(values)
  end
  # println("vals = $(vals)")
  ri.data[keyName] = vals
end

function get_values(ri::RunInfo, keyName)
  ri.data[keyName]
end

function getindex(ri::RunInfo, keyName)
  ri.data[keyName]
end

function setindex!(ri::RunInfo, value, key)
  add_values!(ri, keyName = key, values = value)
end

function print_values(output, keyName )
  open("CMAESoutput.csv", "w") do f
    write(f,"rep,fn,mu,lambda,fit,gen,chr\n")
    for outValues in output
      prm = outValues["parameters"]
      write(f, "$(prm["rep"]),$(prm["fn"]), $(prm["mu"]),$(prm["lambda"]),")
      write(f, "$(out["fit"]),$(out["gen"]),")
      writecsv(f,out["chr"]')
    end
  end
end

function print(runInfo::RunInfo, keyName::Symbol)
  function printvalues(values)
    i = 0
    for value in reverse(values)
      print("$(i):\t"); print(value); print("\n")
      i += 1
    end
  end

  values = runInfo[keyName]
  if (length(values) > 1)
    printvalues(values)
  else
    print(first(values))
  end
end
