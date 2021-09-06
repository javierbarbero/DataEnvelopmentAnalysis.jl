# This file contains function to enable progress meter log 


function progressbarmeter(n::Integer; show_progress = true)
    if show_progress == true 
        Progress(n; enabled = show_progress)
    else
        Progress(n;enabled = false)
    end
end
