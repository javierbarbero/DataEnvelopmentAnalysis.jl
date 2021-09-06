# This file contains function to enable progress meter log 


function progressbarmeter(n::Integer; desc::String = "Computing DEA model", progress::Bool = true)
    if progress == true 
        Progress(n; desc = desc, enabled = progress)
    else
        Progress(n; enabled = false)
    end
end
