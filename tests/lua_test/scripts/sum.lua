---
--- Created by daan.
---

function sum(data)
    local sum = 0

    for i = 1, #data do
        sum = sum + data[i]
    end

    return sum

end